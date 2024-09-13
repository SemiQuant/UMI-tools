import sys
import collections
import re
import os
from builtins import dict
import pysam
import pandas as pd
import numpy as np
import umi_tools
import umi_tools.Utilities as U
import umi_tools.Documentation as Documentation
import umi_tools.network as network
import umi_tools.umi_methods as umi_methods
import umi_tools.sam_methods as sam_methods
import concurrent.futures
import threading
import argparse

# Add the generic docstring text
__doc__ = (__doc__ or "") + Documentation.GENERIC_DOCSTRING_GDC
__doc__ = __doc__ + Documentation.GROUP_DEDUP_GENERIC_OPTIONS

usage = '''
dedup - Deduplicate reads using UMI and mapping coordinates
Usage: umi_tools dedup [OPTIONS] [--stdin=IN_BAM] [--stdout=OUT_BAM]
note: If --stdout is omitted, standard out is output.
To generate a valid BAM file on standard out, please redirect log with --log=LOGFILE or --log2stderr
'''

def detect_bam_features(bamfile, n_entries=1000):
    inbam = pysam.Samfile(bamfile)
    inbam = inbam.fetch(until_eof=True)
    tags = ["NH", "X0", "XT"]
    available_tags = {x: 1 for x in tags}
    for n, read in enumerate(inbam):
        if n > n_entries:
            break
        if read.is_unmapped:
            continue
        else:
            for tag in tags:
                if not read.has_tag(tag):
                    available_tags[tag] = 0
    return available_tags

def aggregateStatsDF(stats_df):
    grouped = stats_df.groupby("UMI")
    agg_dict = {'counts': ['median', len, 'sum']}
    agg_df = grouped.agg(agg_dict)
    agg_df.columns = ['median_counts', 'times_observed', 'total_counts']
    return agg_df

def process_bundle(bundle, key, status, options, processor):
    nOutput = 0
    outreads = []
    stats = {}

    if options.ignore_umi:
        for umi in bundle:
            nOutput += 1
            outreads.append(bundle[umi]["read"])
    else:
        reads, umis, umi_counts = processor(
            bundle=bundle,
            threshold=options.threshold)
        
        if len(reads) == 0:
            return nOutput, outreads, stats

        for read in reads:
            outreads.append(read)
            nOutput += 1

        if options.stats:
            stats['pre'] = {
                'UMI': list(bundle.keys()),
                'counts': [bundle[UMI]['count'] for UMI in bundle]
            }
            stats['post'] = {
                'UMI': umis,
                'counts': umi_counts
            }

    return nOutput, outreads, stats

def main(argv=None):
    if argv is None:
        argv = sys.argv

    # Setup command line parser
    parser = U.OptionParser(version="%prog version: $Id$",
                            usage=usage,
                            description=globals()["__doc__"])

    # Add all options from the original script
    group = U.OptionGroup(parser, "dedup-specific options")
    
    # Add options here (example)
    group.add_option("--output-stats", dest="stats", type="string",
                     default=False, help="Specify location to output stats")
    
    # Add the new threads option
    group.add_option("--threads", dest="threads", type="int",
                     default=1, help="Number of threads to use for processing")

    parser.add_option_group(group)

    # Add common options (-h/--help, ...) and parse command line
    (options, args) = U.Start(parser, argv=argv, add_dedup_count_sam_options=True)

    U.validateSamOptions(options, group=False)

    if options.random_seed:
        np.random.seed(options.random_seed)

    if options.stdin != sys.stdin:
        in_name = options.stdin.name
        options.stdin.close()
    else:
        raise ValueError("Input on standard in not currently supported")

    if options.stdout != sys.stdout:
        if options.no_sort_output:
            out_name = options.stdout.name
        else:
            out_name = U.getTempFilename(dir=options.tmpdir)
            sorted_out_name = options.stdout.name
        options.stdout.close()
    else:
        if options.no_sort_output:
            out_name = "-"
        else:
            out_name = U.getTempFilename(dir=options.tmpdir)
            sorted_out_name = "-"

    if not options.no_sort_output:
        # need to determine the output format for sort
        if options.out_sam:
            sort_format = "sam"
        else:
            sort_format = "bam"

    if options.in_sam:
        in_mode = "r"
    else:
        in_mode = "rb"

    if options.out_sam:
        out_mode = "wh"
    else:
        out_mode = "wb"

    if options.stats and options.ignore_umi:
        raise ValueError("'--output-stats' and '--ignore-umi' options"
                         " cannot be used together")

    infile = pysam.Samfile(in_name, in_mode)
    outfile = pysam.Samfile(out_name, out_mode, template=infile)

    if options.paired:
        outfile = sam_methods.TwoPassPairWriter(infile, outfile)

    nInput, nOutput, input_reads, output_reads = 0, 0, 0, 0

    if options.detection_method:
        bam_features = detect_bam_features(infile.filename)
        if not bam_features[options.detection_method]:
            if sum(bam_features.values()) == 0:
                raise ValueError(
                    "There are no bam tags available to detect multimapping. "
                    "Do not set --multimapping-detection-method")
            else:
                raise ValueError(
                    "The chosen method of detection for multimapping (%s) "
                    "will not work with this bam. Multimapping can be detected"
                    " for this bam using any of the following: %s" % (
                        options.detection_method,
                        ",".join([x for x in bam_features if bam_features[x]])))

    gene_tag = options.gene_tag
    metacontig2contig = None

    if options.chrom:
        inreads = infile.fetch(reference=options.chrom)
    else:
        if options.per_contig and options.gene_transcript_map:
            metacontig2contig = sam_methods.getMetaContig2contig(
                infile, options.gene_transcript_map)
            metatag = "MC"
            inreads = sam_methods.metafetcher(infile, metacontig2contig, metatag)
            gene_tag = metatag
        else:
            inreads = infile.fetch()

    # set up ReadCluster functor with methods specific to
    # specified options.method
    processor = network.ReadDeduplicator(options)

    bundle_iterator = sam_methods.get_bundles(
        options, metacontig_contig=metacontig2contig)

    if options.stats:
        # set up arrays to hold stats data
        stats_pre_df_dict = {"UMI": [], "counts": []}
        stats_post_df_dict = {"UMI": [], "counts": []}
        pre_cluster_stats = []
        post_cluster_stats = []
        pre_cluster_stats_null = []
        post_cluster_stats_null = []
        topology_counts = collections.Counter()
        node_counts = collections.Counter()

    read_gn = umi_methods.random_read_generator(
        infile.filename,
        chrom=options.chrom,
        barcode_getter=bundle_iterator.barcode_getter)

    # Create a thread-safe list for output
    output_list = []
    output_lock = threading.Lock()

    # Create a ThreadPoolExecutor
    with concurrent.futures.ThreadPoolExecutor(max_workers=options.threads) as executor:
        futures = []

        for bundle, key, status in bundle_iterator(inreads):
            nInput += sum([bundle[umi]["count"] for umi in bundle])

            future = executor.submit(process_bundle, bundle, key, status, options, processor)
            futures.append(future)

            # Process completed futures
            for completed in concurrent.futures.as_completed(futures):
                nOutput_bundle, outreads_bundle, stats_bundle = completed.result()
                nOutput += nOutput_bundle

                with output_lock:
                    output_list.extend(outreads_bundle)

                if options.stats:
                    # Update stats
                    if 'pre' in stats_bundle:
                        stats_pre_df_dict['UMI'].extend(stats_bundle['pre']['UMI'])
                        stats_pre_df_dict['counts'].extend(stats_bundle['pre']['counts'])
                    if 'post' in stats_bundle:
                        stats_post_df_dict['UMI'].extend(stats_bundle['post']['UMI'])
                        stats_post_df_dict['counts'].extend(stats_bundle['post']['counts'])

                # Write processed reads to file
                while len(output_list) >= 100000:
                    with output_lock:
                        for _ in range(100000):
                            outfile.write(output_list.pop(0))
                    output_reads += 100000
                    U.info("Written out %i reads" % output_reads)

                while nInput >= input_reads + 1000000:
                    input_reads += 1000000
                    U.info("Parsed %i input reads" % input_reads)

        # Write remaining reads
        with output_lock:
            for read in output_list:
                outfile.write(read)
            nOutput += len(output_list)

    outfile.close()

    if not options.no_sort_output:
        # Sort the output
        pysam.sort("-o", sorted_out_name, "-O", sort_format, "--no-PG", out_name)
        os.unlink(out_name)  # Delete the tempfile

    if options.stats:
        # generate the stats dataframe
        stats_pre_df = pd.DataFrame(stats_pre_df_dict)
        stats_post_df = pd.DataFrame(stats_post_df_dict)

        # tally the counts per umi per position
        pre_counts = collections.Counter(stats_pre_df["counts"])
        post_counts = collections.Counter(stats_post_df["counts"])

        counts_index = list(set(pre_counts.keys()).union(set(post_counts.keys())))
        counts_index.sort()

        with U.openFile(options.stats + "_per_umi_per_position.tsv", "w") as outf:
            outf.write("counts\tinstances_pre\tinstances_post\n")
            for count in counts_index:
                values = (count, pre_counts[count], post_counts[count])
                outf.write("\t".join(map(str, values)) + "\n")

        # aggregate stats pre/post per UMI
        agg_pre_df = aggregateStatsDF(stats_pre_df)
        agg_post_df = aggregateStatsDF(stats_post_df)

        agg_df = pd.merge(agg_pre_df, agg_post_df, how='left',
                          left_index=True, right_index=True,
                          sort=True, suffixes=["_pre", "_post"])

        # TS - if count value not observed either pre/post-dedup,
        # merge will leave an empty cell and the column will be cast as a float
        # see http://pandas.pydata.org/pandas-docs/dev/missing_data.html
        # --> Missing data casting rules and indexing
        # so, back fill with zeros and convert back to int
        agg_df = agg_df.fillna(0).astype(int)

        agg_df.index = [x.decode() for x in agg_df.index]
        agg_df.index.name = 'UMI'
        agg_df.to_csv(options.stats + "_per_umi.tsv", sep="\t")

        # bin distances into integer bins
        max_ed = int(max(map(max, [pre_cluster_stats, post_cluster_stats,
                                   pre_cluster_stats_null, post_cluster_stats_null])))

        cluster_bins = range(-1, int(max_ed) + 2)

        def bin_clusters(cluster_list, bins=cluster_bins):
            ''' take list of floats and return bins'''
            return np.digitize(cluster_list, bins, right=True)

        def tallyCounts(binned_cluster, max_edit_distance):
            ''' tally counts per bin '''
            return np.bincount(binned_cluster,
                               minlength=max_edit_distance + 3)

        pre_cluster_binned = bin_clusters(pre_cluster_stats)
        post_cluster_binned = bin_clusters(post_cluster_stats)
        pre_cluster_null_binned = bin_clusters(pre_cluster_stats_null)
        post_cluster_null_binned = bin_clusters(post_cluster_stats_null)

        edit_distance_df = pd.DataFrame({
            "unique": tallyCounts(pre_cluster_binned, max_ed),
            "unique_null": tallyCounts(pre_cluster_null_binned, max_ed),
            options.method: tallyCounts(post_cluster_binned, max_ed),
            "%s_null" % options.method: tallyCounts(post_cluster_null_binned, max_ed),
            "edit_distance": cluster_bins
        }, columns=["unique", "unique_null", options.method,
                    "%s_null" % options.method, "edit_distance"])

        edit_distance_df['edit_distance'] = edit_distance_df['edit_distance'].astype(str)
        edit_distance_df.loc[0, 'edit_distance'] = "Single_UMI"

        edit_distance_df.to_csv(options.stats + "_edit_distance.tsv",
                                index=False, sep="\t")

    # write footer and output benchmark information.
    U.info("Reads: %s" % ", ".join(
        ["%s: %s" % (x[0], x[1]) for x in
         bundle_iterator.read_events.most_common()]))

    U.info("Number of reads out: %i" % nOutput)

    if not options.ignore_umi:
        U.info("Total number of positions deduplicated: %i" %
               processor.UMIClusterer.positions)

        if processor.UMIClusterer.positions > 0:
            U.info("Mean number of unique UMIs per position: %.2f" %
                   (float(processor.UMIClusterer.total_umis_per_position) /
                    processor.UMIClusterer.positions))
            U.info("Max. number of unique UMIs per position: %i" %
                    processor.UMIClusterer.max_umis_per_position)
    else:
        U.warn("The BAM did not contain any valid "
               "reads/read pairs for deduplication")

    if options.filter_umi:
        U.info("%i UMIs were in a group where the top UMI was not a "
               "whitelist UMI and were therefore discarded" %
               processor.umi_whitelist_counts["Non-whitelist UMI"])

    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
