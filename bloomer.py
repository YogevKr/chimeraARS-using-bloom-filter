import argparse
import math
from itertools import islice

from Bio import SeqIO
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


def count_hits(k_bloom_tuple, target_genome):
    k, bloom = k_bloom_tuple
    return sum(map(lambda x: x in bloom, window(target_genome, n=k)))


def chimera_ars_score(
    host_genome: str,
    target_genome: str,
    max_k: int = 20,
    error_rate: float = 0.0001,
    num_processors: int = 0,
) -> float:
    from multiprocessing import Pool, cpu_count
    import functools

    pool = Pool(num_processors or cpu_count())

    populate_runner = functools.partial(
        get_populted_bloom_filter, host_genome=host_genome, error_rate=error_rate
    )
    bloom_models = pool.map(populate_runner, range(1, max_k + 1))

    hit_runner = functools.partial(count_hits, target_genome=target_genome)
    total_hits_list = pool.map(
        hit_runner, ((k, bloom) for k, bloom in enumerate(bloom_models, start=1))
    )

    print(total_hits_list)

    return sum(total_hits_list) / len(target_genome)


def main(args):
    host_genome: str = tuple(islice(SeqIO.parse(args.host_genome, "fasta"), 1))[0]
    target_genome: str = tuple(islice(SeqIO.parse(args.target_genome, "fasta"), 1))[0]

    score = chimera_ars_score(
        host_genome=host_genome.seq,
        target_genome=target_genome.seq,
        max_k=args.max_k,
        error_rate=args.error_rate / args.max_k,
        num_processors=0,
    )
    print(score)


if __name__ == "__main__":
    # Parse argument
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--host_genome", type=str, required=True, help="Path of host genome FASTA file"
    )
    parser.add_argument(
        "--target_genome", type=str, required=True, help="Path of target genome FASTA file"
    )
    parser.add_argument("--max_k", type=int, default=20, help="Longest sequence")
    parser.add_argument("--error_rate", type=int, default=0.01, help="Error rate")
    parser.add_argument(
        "--num_processors", type=int, default=0, help="Number of processors"
    )

    args = parser.parse_args()
    main(args)
