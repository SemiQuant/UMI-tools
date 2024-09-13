'''
===========================================================
dedup - Deduplicate reads using UMI and mapping coordinates
===========================================================

*Deduplicate reads based on the mapping co-ordinate and the UMI attached to the read*

The identification of duplicate reads is performed in an error-aware
manner by building networks of related UMIs (see
``--method``). ``dedup`` can also handle cell barcoded input (see
``--per-cell``).

Usage::

    umi_tools dedup --stdin=INFILE --log=LOGFILE [OPTIONS] > OUTFILE

Selecting the representative read
---------------------------------
For every group of duplicate reads, a single representative read is
retained.The following criteria are applied to select the read that
will be retained from a group of duplicated reads:

1. The read with the lowest number of mapping coordinates (see
``--multimapping-detection-method`` option)

2. The read with the highest mapping quality. Note that this is not
the read sequencing quality and that if two reads have the same
mapping quality then one will be picked at random regardless of the
read quality.

Otherwise a read is chosen at random.


Dedup-specific options
----------------------
"""""""""""""""""""""""""""
``--output-stats=[PREFIX]``
"""""""""""""""""""""""""""
One can use the edit distance between UMIs at the same position as an
quality control for the deduplication process by comparing with
a null expectation of random sampling. For the random sampling, the
observed frequency of UMIs is used to more reasonably model the null
expectation.

Use this option to generate a stats outfile called:

[PREFIX]_edit_distance.tsv
  Reports the (binned) average edit distance between the UMIs at each
  position. Positions with a single UMI are reported seperately.  The
  edit distances are reported pre- and post-deduplication alongside
  the null expectation from random sampling of UMIs from the UMIs
  observed across all positions. Note that separate null
  distributions are reported since the null depends on the observed
  frequency of each UMI which is different pre- and
  post-deduplication. The post-duplication values should be closer to
  their respective null than the pre-deduplication vs null comparison

In addition, this option will trigger reporting of further summary
statistics for the UMIs which may be informative for selecting the
optimal deduplication method or debugging.

Each unique UMI sequence may be observed [0-many] times at multiple
positions in the BAM. The following files report the distribution for
the frequencies of each UMI.

[PREFIX]_per_umi_per_position.tsv
  The `_per_umi_per_position.tsv` file simply tabulates the
  counts for unique combinations of UMI and position. E.g if prior to
  deduplication, we have two positions in the BAM (POSa, POSb), at
  POSa we have observed 2*UMIa, 1*UMIb and at POSb: 1*UMIc, 3*UMId,
  then the stats file is populated thus:

  ====== =============
  counts instances_pre
  ------ -------------
  1      2
  2      1
  3      1
  ====== =============


  If post deduplication, UMIb is grouped with UMIa such that POSa:
  3*UMIa, then the `instances_post` column is populated thus:

  ====== ============= ==============
  counts instances_pre instances_post
  ------ ------------- --------------
  1      2             1
  2      1             0
  3      1             2
  ====== ============= ==============

[PREFIX]_per_umi.tsv
  The `_per_umi.tsv` table provides UMI-level summary
  statistics. Keeping in mind that each unique UMI sequence can be
  observed at [0-many] times across multiple positions in the BAM,

  :times_observed: How many positions the UMI was observed at
  :total_counts: The total number of times the UMI was observed across all positions
  :median_counts: The median for the distribution of how often the UMI was observed at                  each position (excluding zeros)

  Hence, whenever times_observed=1, total_counts==median_counts.

'''

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

# Add the generic docstring text
__doc__ = __doc__ + Documentation.GENERIC_DOCSTRING_GDC
__doc__ = __doc__ + Documentation.GROUP_DEDUP_GENERIC_OPTIONS

usage = '''
dedup - Deduplicate reads using UMI and mapping coordinates
Usage: umi_tools dedup [OPTIONS] [--stdin=IN_BAM] [--stdout=OUT_BAM]
note: If --stdout is omitted, standard out is output.
To generate a valid BAM file on standard out, please redirect log with --log=LOGFILE or --log2stderr
'''

def detect_bam_features(bamfile, n_entries=1000):
    # ... (unchanged)

def aggregateStatsDF(stats_df):
    # ... (unchanged)

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

    # ... (Add all options as in the original script)

    (options, args) = U.Start(parser, argv=argv, add_dedup_count_sam_options=True)

    # ... (Validation and setup code remains the same)

    infile = pysam.Samfile(in_name, in_mode)
    outfile = pysam.Samfile(out_name, out_mode, template=infile)

    if options.paired:
        outfile = sam_methods.TwoPassPairWriter(infile, outfile)

    nInput, nOutput, input_reads, output_reads = 0, 0, 0, 0

    # ... (Setup code remains the same)

    processor = network.ReadDeduplicator(options)
    bundle_iterator = sam_methods.get_bundles(
        options, metacontig_contig=metacontig2contig)

    if options.stats:
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
    with concurrent.futures.ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
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
        # Generate stats (this part remains mostly unchanged)
        # ... (Stats generation code)

    # Write footer and output benchmark information
    U.info("Reads: %s" % ", ".join(["%s: %s" % (x[0], x[1]) for x in bundle_iterator.read_events.most_common()]))
    U.info("Number of reads out: %i" % nOutput)

    if not options.ignore_umi:
        U.info("Total number of positions deduplicated: %i" % processor.UMIClusterer.positions)
        if processor.UMIClusterer.positions > 0:
            U.info("Mean number of unique UMIs per position: %.2f" % (float(processor.UMIClusterer.total_umis_per_position) / processor.UMIClusterer.positions))
            U.info("Max. number of unique UMIs per position: %i" % processor.UMIClusterer.max_umis_per_position)
    else:
        U.warn("The BAM did not contain any valid reads/read pairs for deduplication")

    if options.filter_umi:
        U.info("%i UMIs were in a group where the top UMI was not a whitelist UMI and were therefore discarded" % processor.umi_whitelist_counts["Non-whitelist UMI"])

    U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
