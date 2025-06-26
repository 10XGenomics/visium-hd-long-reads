#!/usr/bin/env python

"""Process long-reads data to pseudo-short read fastqs compatible to run Space Ranger."""

## Example Usage

# long_reads_to_10x_paired.py --bam <input_long-reads.bam> --sample_name <sample_name> [--optional flags]

import argparse
import gzip
import json
import multiprocessing as mp
import queue
import threading
import time
from collections import Counter, namedtuple
from collections.abc import Iterator
from functools import partial
from queue import Queue
from typing import Any

import edlib
from pysam import AlignmentFile, FastxFile

start_time = time.perf_counter()

complement_trans = str.maketrans("ACGTacgt", "TGCAtgca")

# Visium HD 3' Adapter Sequences
adapters = {
    "adapter1_f": "CTACACGACGCTCTTCCGATCT",
    "adapter1_r": "AGATCGGAAGAGCGTCGTGTAG",
    "adapter2_f": "ATGTACTCTGCGTTGATACCACTGCTT",
    "adapter2_r": "AAGCAGTGGTATCAACGCAGAGTACAT",
}

valid_adapter_pairs = {
    ("adapter1_f", "adapter2_f"): "f",
    ("adapter2_r", "adapter1_r"): "r",
}

Read = namedtuple("Read", ["name", "sequence", "quality"])

VISIUM_HD_R1_LENGTH = 43


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments.

    Returns:
        argparse.Namespace: parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Generate short-read R1/R2 fastqs from long-read BAM/fastqs."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--fastq", help="Path to the input long-read .fastq file.")
    group.add_argument("--bam", help="Path to the input long read .bam file.")
    parser.add_argument("--threads", default=12, type=int, help="Threads to use.")
    parser.add_argument(
        "--compress",
        default=True,
        type=bool,
        help="Compression of the R1/R2 .fastq outputs.",
    )
    parser.add_argument(
        "--chunk_size", default=2500, type=int, help="Chunksize for multiprocessing."
    )
    parser.add_argument(
        "--sample_name", default="SAMPLE", help="Sample name for output file prefix."
    )
    parser.add_argument(
        "--r1_size",
        default=43,
        type=int,
        help="Number of bases after the adapter to include in R1 fastq file.",
    )
    parser.add_argument(
        "--r2_size",
        default=200,
        type=int,
        help="Number of bases after the adapter to include in R2 fastq file.",
    )
    parser.add_argument(
        "--min_id", default=0.8, type=float, help="Minimum identity to call an adapter."
    )
    return parser.parse_args()


def read_chunks(
    input_file: str, input_format: str, chunk_size: int = 50000
) -> Iterator[list[Read]]:
    """Read the input long-reads data in chunks."""
    chunk = []

    if input_format == "bam":
        # Pysam AlignmentFile objects are not pickleable, so we need to convert the data
        # into a pickleable format.
        with AlignmentFile(input_file, check_sq=False, threads=2) as fh:
            for record in fh.fetch(until_eof=True):
                read = Read(
                    record.query_name,
                    record.query_sequence,
                    "".join(
                        chr(q + 33) for q in record.query_qualities  # type: ignore[possibly-none]
                    ),  ## Convert quality scores to ASCII
                )
                chunk.append(read)
                if len(chunk) == chunk_size:
                    yield chunk
                    chunk = []
            if chunk:
                yield chunk
    else:
        # pysam FastxFile objects appear to be pickleable, so we can use them with multiprocessing
        with FastxFile(input_file) as fh:
            for record in fh:
                chunk.append(record)
                if len(chunk) == chunk_size:
                    yield chunk
                    chunk = []
            if chunk:
                yield chunk


def get_subreads(
    read: Read, r1_size: int, r2_size: int, min_id: float
) -> tuple[list[str], list[str], list[str], list[str]]:
    """Split long-reads into R1/R2 sub-reads based on presence of compatible adapter pairs."""
    r1 = []
    r2 = []
    hit_summary = []
    config_summary = []
    edlib_hits = []
    valid_edlib_hits = []

    for adapter, adapter_seq in adapters.items():
        max_ed = int(len(adapter_seq) * (1 - min_id))
        result = edlib.align(
            adapter_seq, read.sequence, mode="HW", task="locations", k=max_ed
        )

        if result["editDistance"] != -1 and result["editDistance"] <= max_ed:
            for loc in result["locations"]:
                edlib_hits.append([adapter, loc])

    # Sort the hits by start position
    edlib_hits.sort(key=lambda x: x[1][0])

    config_summary.append("-".join([x[0] for x in edlib_hits]))

    n_valid_segments = 0
    for i in range(len(edlib_hits) - 1):
        # Get adapters that would be first in a valid pair, depending on strand
        if edlib_hits[i][0] not in ["adapter1_f", "adapter2_r"]:
            continue

        a1 = edlib_hits[i]
        a2 = edlib_hits[i + 1]
        subread_orientation = valid_adapter_pairs.get((a1[0], a2[0]))
        if subread_orientation is None:
            continue

        a1_end = a1[1][1] + 1
        a2_start = a2[1][0]
        if abs(a1_end - a2_start) < VISIUM_HD_R1_LENGTH:
            continue

        valid_edlib_hits.append([a1, a2, subread_orientation])
        n_valid_segments += 1

    for i, (a1, a2, subread_orientation) in enumerate(valid_edlib_hits):
        subread_id = (
            f"{read.name}" if len(valid_edlib_hits) == 1 else f"{read.name}_{i}"
        )

        a1_end = a1[1][1] + 1
        a2_start = a2[1][0]
        subread_seq = read.sequence[a1_end:a2_start]
        subread_qual = read.quality[a1_end:a2_start]

        if subread_orientation == "r":
            # For Visium HD 3', reads with F adapters will be in the mRNA
            # antisense orientation and need to be reverse-complemented so
            # that output reads are in mRNA sense orientation.
            subread_seq = subread_seq[::-1].translate(complement_trans)
            subread_qual = subread_qual[::-1]

        hit_summary.append(
            f"{subread_id}\t{subread_orientation}\t"
            f"{a1[0]}\t{a1_end}\t{a2[0]}\t{a2_start}\n"
        )

        r1.append(
            f"@{subread_id} 1:N:0:0\n{subread_seq[:r1_size]}\n+\n{subread_qual[:r1_size]}\n"
        )
        # r2 reads are expected in opposite orientation to r1
        r2.append(
            f"@{subread_id} 4:N:0:0\n{subread_seq[-r2_size:][::-1].translate(complement_trans)}\n+\n{subread_qual[-r2_size:][::-1]}\n"
        )

    return r1, r2, hit_summary, config_summary


def process_chunk(
    chunk: list[Read],
    r1_size: int,
    r2_size: int,
    min_id: float,
    r1_q: Queue[Any],
    r2_q: Queue[Any],
    summary_q: Queue[Any],
) -> Counter:
    """Process the input chunks and add to queues. Returns a Counter object for aggregation."""
    local_config_counter = Counter()
    r1_chunk = []
    r2_chunk = []
    summary_chunk = []

    for read in chunk:
        (
            r1_res,
            r2_res,
            summ,
            configs,
        ) = get_subreads(read, r1_size, r2_size, min_id)
        if r1_res:
            r1_chunk.extend(r1_res)
            r2_chunk.extend(r2_res)
            summary_chunk.extend(summ)
        if configs:
            local_config_counter.update(configs)
    if r1_chunk:
        r1_q.put(r1_chunk)
        r2_q.put(r2_chunk)
        summary_q.put(summary_chunk)

    return local_config_counter


def writer_thread(write_queue: queue.Queue, r_path: str, compress: bool = True) -> None:
    """Write FASTQ entries from queue to files."""
    if compress:
        f = gzip.open(r_path, "wt", compresslevel=1)
    else:
        f = open(r_path, "w")

    with f as r_fh:
        while True:
            item = write_queue.get()
            if item is None:
                break
            r_fh.writelines(item)


def main(args: argparse.Namespace) -> None:
    """Extract R1 and R2 reads, splitting chimeric reads if necessary."""
    input_fmt = "bam" if args.bam else "fastq"

    config_counter = Counter()

    pool = mp.Pool(args.threads)
    manager = mp.Manager()
    r1_write_queue = manager.Queue(maxsize=2 * args.threads)
    r2_write_queue = manager.Queue(maxsize=2 * args.threads)
    summary_write_queue = manager.Queue(maxsize=2 * args.threads)

    worker_func = partial(
        process_chunk,
        r1_size=args.r1_size,
        r2_size=args.r2_size,
        min_id=args.min_id,
        r1_q=r1_write_queue,
        r2_q=r2_write_queue,
        summary_q=summary_write_queue,
    )

    ext = "fastq.gz" if args.compress else "fastq"
    # Start writer threads
    r1_fname = f"{args.sample_name}_S1_L001_R1_001.{ext}"
    r2_fname = f"{args.sample_name}_S1_L001_R2_001.{ext}"
    r1_writer = threading.Thread(
        target=writer_thread, args=(r1_write_queue, r1_fname, args.compress)
    )
    r2_writer = threading.Thread(
        target=writer_thread, args=(r2_write_queue, r2_fname, args.compress)
    )
    summary_writer = threading.Thread(
        target=writer_thread,
        args=(summary_write_queue, f"{args.sample_name}_summary.tsv.gz", True),
    )
    r1_writer.start()
    r2_writer.start()
    summary_writer.start()

    input_file = args.bam if args.bam else args.fastq

    for config_counter_chunk in pool.imap_unordered(
        worker_func,
        read_chunks(input_file, input_fmt, args.chunk_size),
    ):
        config_counter.update(config_counter_chunk)

    # Signal the writer threads to exit
    r2_write_queue.put(None)
    r1_write_queue.put(None)
    summary_write_queue.put(None)
    r1_writer.join()
    r2_writer.join()
    summary_writer.join()

    # Write out summary counts of the adapter configurations to configs.json
    with open(f"{args.sample_name}.configs.json", "w") as f:
        json.dump(config_counter, f, indent=4)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)

finish_time = time.perf_counter()
print(f"Finished in {round(finish_time - start_time, 2)} second(s)")
