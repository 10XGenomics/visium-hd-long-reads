#!/usr/bin/env python

"""Add 10x Space Ranger BC/UMI and spatial tags to raw long-read .bams."""

## Usage

# add_10x_bam_tags.py --sr_bam <possorted_genome_bam.bam> --lr_bam <long-reads.bam> [--optional flags]

## Libraries
import argparse
import csv
import gzip
import os
import queue
import threading
import time
from array import array
from collections import defaultdict
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor, as_completed

import pysam

TagValue = str | int | float | array | None

start_time = time.perf_counter()


def cli_args() -> argparse.Namespace:
    """Specify command-line arguments for inputs and outputs."""
    parser = argparse.ArgumentParser(
        description="Tag raw long-reads .bam with 10x Space Ranger Visium HD BC/UMI and spatial coordinates"
    )
    parser.add_argument(
        "--sr_bam",
        required=True,
        type=str,
        help="Path to the Space Ranger output possorted_genome_bam.bam file.",
    )
    parser.add_argument(
        "--lr_bam", required=True, type=str, help="Path to the raw long-reads BAM file"
    )
    parser.add_argument(
        "--chunk_size",
        type=int,
        default=20000,
        help="Optional chunk size (default: 20000).",
    )
    parser.add_argument(
        "--out_dir",
        required=False,
        type=str,
        default=os.getcwd(),
        help="Specify output path (default: current working directory)",
    )
    parser.add_argument(
        "--overwrite",
        required=False,
        action="store_true",
        help="Boolean to optionally overwrite the raw long-reads BAM as output (default: False)",
    )

    args = parser.parse_args()

    print(f"Space Ranger BAM file: {args.sr_bam}")
    print(f"Long-Reads BAM file: {args.lr_bam}")
    print(f"Chunk size: {args.chunk_size}")
    print(f"Output directory: {args.out_dir}")
    print(f"Overwrite Long-Reads BAM: {args.overwrite}")

    return args


if __name__ == "__cli_args__":
    cli_args()


def chunk_input_bam(
    sr_bam: str, chunk_size: int, output_csv_gz: str, thread_count: int = 4
) -> Iterator[defaultdict[str | None, dict[str, TagValue]]]:
    """Read the processed Space Ranger .bam in chunks and get read tags.

    Args:
        sr_bam (bam): Space Ranger output; possorted_genome_bam.bam
        chunk_size (int): size to chunk .bam
        output_csv_gz (file): output csv filepath
        thread_count (int): thread count to use
    """
    chunk_read_tags: defaultdict[str | None, dict[str, TagValue]] = defaultdict(
        lambda: {"CR": None, "CB": None, "UR": None, "UB": None}
    )

    with (
        pysam.AlignmentFile(
            sr_bam, "rb", check_sq=False, threads=thread_count
        ) as processed_bam,
        gzip.open(output_csv_gz, "wt", newline="") as csv_gz,
    ):

        csv_writer = csv.writer(csv_gz)
        csv_writer.writerow(
            [
                "read_name",
                "uncorrected_barcode",
                "corrected_barcode",
                "uncorrected_umi",
                "corrected_umi",
            ]
        )

        for read in processed_bam.fetch(until_eof=True):

            cr_tag = read.get_tag("CR") if read.has_tag("CR") else None
            cb_tag = read.get_tag("CB") if read.has_tag("CB") else None
            ur_tag = read.get_tag("UR") if read.has_tag("UR") else None
            ub_tag = read.get_tag("UB") if read.has_tag("UB") else None

            ## Add read/BC/UMI info to chunk dict
            chunk_read_tags[read.query_name].update(
                {"CR": cr_tag, "CB": cb_tag, "UR": ur_tag, "UB": ub_tag}
            )
            ## Also write to .csv.gz
            csv_writer.writerow([read.query_name, cr_tag, cb_tag, ur_tag, ub_tag])

            if len(chunk_read_tags) == chunk_size:
                yield chunk_read_tags
                chunk_read_tags = defaultdict(
                    lambda: {"CR": None, "CB": None, "UR": None, "UB": None}
                )

        if chunk_read_tags:
            yield chunk_read_tags


def writer_thread_func(bam: pysam.AlignmentFile, q: queue.Queue) -> None:
    """Writer thread function to write reads from queue to BAM file."""
    while True:
        reads = q.get()
        if reads is None:  # A 'None' object is the sentinel to stop the thread.
            break
        for read in reads:
            bam.write(read)
        q.task_done()


def process_read_chunk(
    reads: list[pysam.AlignedSegment],
    processed_tags: dict[str | None, dict[str, TagValue]],
    q: queue.Queue,
) -> int:
    """Processes a chunk of reads, adds tags, and puts the chunk on the writer queue."""
    tagged_reads = 0
    processed_reads = []
    for read in reads:
        if read.query_name in processed_tags:
            tagged_reads += 1
            for tag, value in processed_tags[read.query_name].items():
                if value is not None:
                    read.set_tag(tag, value, value_type="Z")
        processed_reads.append(read)
    q.put(processed_reads)
    return tagged_reads


def update_bam_tags(
    chunk_input_bam_generator: Iterator[defaultdict[str | None, dict[str, TagValue]]],
    lr_bam: str,
    tagged_bam: str,
    overwrite_bam: bool,
    chunk_size: int,
    thread_count: int = 4,
) -> None:
    """Consume chunk_input_bam generator, iterate over the raw long-reads bam to update bam tags.

    Args:
        chunk_input_bam_generator (dict): generator of dictionaries with read_name: {tags}
        lr_bam (bam): raw long-reads input bam to update tags
        tagged_bam (bam): output .bam with BC/UMI tags
        overwrite_bam (bool): if true, overwrite lr_bam input bam
        chunk_size (int): number of records to chunk bam
        thread_count (int): thread count to use
    """
    # Create a dictionary to store read names and new tag values
    processed_bam_tags = {
        read_name: tags
        for chunk in chunk_input_bam_generator
        for read_name, tags in chunk.items()
    }
    print(f"Reads in Space Ranger BAM: {len(processed_bam_tags)}")

    # Setup up queue and writer thread
    writer_queue = queue.Queue(maxsize=thread_count * 2)

    with (
        pysam.AlignmentFile(
            lr_bam, "rb", check_sq=False, threads=thread_count
        ) as raw_bam,
        pysam.AlignmentFile(tagged_bam, "wb", header=raw_bam.header) as out_bam,
        ThreadPoolExecutor(max_workers=thread_count) as executor,
    ):

        # Start the writer thread
        writer_thread = threading.Thread(
            target=writer_thread_func, args=(out_bam, writer_queue)
        )
        writer_thread.daemon = (
            True  # Allows main thread to exit if this thread is blocked
        )
        writer_thread.start()

        # This set holds the 'Future' objects for tasks that are currently running or queued to bound the work queue.
        futures = set()
        max_futures_in_flight = thread_count * 2

        input_reads = 0
        output_reads_tagged = 0

        bam_iterator = raw_bam.fetch(until_eof=True)
        finished_reading = False

        # This loop continues as long as we are reading the input BAM or there are tasks in flight.
        while not finished_reading or futures:
            # Stage 1: Submit chunks depending on available futures and reading status.
            while len(futures) < max_futures_in_flight and not finished_reading:
                chunk = []
                try:
                    for _ in range(chunk_size):
                        chunk.append(next(bam_iterator))
                except StopIteration:
                    # This indicates we've reached the end of the input BAM file.
                    finished_reading = True

                if chunk:
                    input_reads += len(chunk)
                    future = executor.submit(
                        process_read_chunk, chunk, processed_bam_tags, writer_queue
                    )
                    futures.add(future)
                # If chunk is empty and reading is finished, exit loop.

            # Stage 2: Process completed tasks
            for future in as_completed(futures):
                output_reads_tagged += future.result()
                # Remove the completed future from our tracking set.
                futures.remove(future)
                # Break after processing completed future to submit a new chunk to queue.
                break
        # All chunks have been processed. Signal the writer thread to stop and wait for it.
        writer_queue.put(None)
        writer_thread.join()

    # Optionally replace the old raw.bam file with the updated one
    if overwrite_bam:
        os.replace(tagged_bam, lr_bam)

    print(f"Long reads processed: {input_reads}")
    print(f"Long reads tagged: {output_reads_tagged}")
    if output_reads_tagged < input_reads:
        print(
            "(Long reads that were split into subreads will not be tagged but are present in the CSV output)"
        )


def main() -> None:
    """Add Space Ranger BAM tags to long-read BAM files."""
    args = cli_args()

    ## specify output filenames
    output_csv_gz = os.path.join(
        args.out_dir,
        f'{os.path.basename(args.lr_bam).rsplit(".bam", 1)[0]}.spatial_barcodes.csv.gz',
    )
    tagged_bam = os.path.join(
        args.out_dir, f'{os.path.basename(args.lr_bam).rsplit(".bam", 1)[0]}.tagged.bam'
    )

    chunk_generator = chunk_input_bam(args.sr_bam, args.chunk_size, output_csv_gz)
    update_bam_tags(
        chunk_generator, args.lr_bam, tagged_bam, args.overwrite, args.chunk_size
    )


if __name__ == "__main__":
    main()

finish_time = time.perf_counter()
print(f"Finished in {round(finish_time - start_time, 2)} second(s)")
