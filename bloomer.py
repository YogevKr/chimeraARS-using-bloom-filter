import argparse
import math
from itertools import islice

from Bio import SeqIO, Seq
from bloom_filter import BloomFilter
from tqdm import tqdm


def window(seq, n=2):
    """
    Returns a sliding window (of width n) over data from the iterable
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   
    """
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result


def num_of_bits(ideal_num_elements_n, error_rate):
    numerator = -1 * ideal_num_elements_n * math.log(error_rate)
    denominator = math.log(2) ** 2
    return numerator / denominator


def get_populted_bloom_filter(k, host_genome, error_rate):
    max_elements = 4 ** k if 4 ** k < len(host_genome) else len(host_genome)

    bloom = BloomFilter(max_elements=max_elements, error_rate=error_rate)
    print(f"k={k}, max_elements={max_elements}, num_bits_m={bloom.num_bits_m}")

    for s in window(host_genome, n=k):
        bloom.add("".join(s))

    return bloom


def count_hits(bloom, k, target_genome):
    return sum(map(lambda x: x in bloom, window(target_genome, n=k)))


def run_one_window(k, *, host_genome, target_genome, error_rate):
    bloom = get_populted_bloom_filter(k, host_genome, error_rate)

    return count_hits(bloom, k, target_genome)


def chimera_ars_score(
    host_genome: Seq.Seq,
    target_genome: Seq.Seq,
    max_k: int = 20,
    error_rate: float = 0.0001,
    num_processors: int = 0,
) -> float:
    from multiprocessing import Pool, cpu_count
    import functools

    pool = Pool(num_processors or cpu_count())

    runner = functools.partial(
        run_one_window,
        host_genome=host_genome,
        target_genome=target_genome,
        error_rate=error_rate,
    )
    total_hits_list = pool.map(runner, range(1, max_k + 1))

    print(total_hits_list)

    return sum(total_hits_list) / len(target_genome)


def main(args):
    host_genome_seq: Seq.Seq = sum(
        [g.seq for g in SeqIO.parse(args.host_genome, "fasta")], Seq.Seq("")
    ).upper()
    target_genome_seq: Seq.Seq = tuple(
        islice(SeqIO.parse(args.target_genome, "fasta"), 1)
    )[0].seq.upper()

    score = chimera_ars_score(
        host_genome=host_genome_seq,
        target_genome=target_genome_seq,
        max_k=args.max_k,
        error_rate=args.error_rate / args.max_k,
        num_processors=0,
    )

    false_positive_rate = 1 - ((score - args.error_rate) / score)

    print(f"Host genome length: {len(host_genome_seq)}")
    print(f"Target genome length: {len(target_genome_seq)}")
    print(
        f"chimeraARS score: {score} ({false_positive_rate * 100:0.2f}% false positive)"
    )


if __name__ == "__main__":
    # Parse argument
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--host_genome", type=str, required=True, help="Path of host genome FASTA file"
    )
    parser.add_argument(
        "--target_genome",
        type=str,
        required=True,
        help="Path of target genome FASTA file (Single sequence)",
    )
    parser.add_argument("--max_k", type=int, default=20, help="Longest sequence")
    parser.add_argument("--error_rate", type=float, default=0.01, help="Error rate")
    parser.add_argument(
        "--num_processors", type=int, default=0, help="Number of processors"
    )

    args = parser.parse_args()
    main(args)
