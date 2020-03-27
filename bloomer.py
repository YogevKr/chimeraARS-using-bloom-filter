from bloom_filter import BloomFilter
from itertools import islice
from tqdm import tqdm
from Bio import SeqIO
import math


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
    max_elements = 4**k if 4**k < len(host_genome) else len(host_genome)

    bloom = BloomFilter(max_elements=max_elements, error_rate=error_rate)
    print(f"k={k}, max_elements={max_elements}, num_bits_m={bloom.num_bits_m}")

    for s in window(host_genome, n=k):
        bloom.add("".join(s))

    return bloom


def count_hits(k_bloom_tuple, target_genome):
    k, bloom = k_bloom_tuple
    return sum(map(lambda x: x in bloom, window(target_genome, n=k)))


def chimera_ars_score(host_genome, target_genome, max_k=16, error_rate=0.0001, num_processors=0):
    from multiprocessing import Pool, cpu_count
    import functools

    pool = Pool(num_processors or cpu_count())

    populate_runner = functools.partial(get_populted_bloom_filter, host_genome=host_genome, error_rate=error_rate)
    bloom_models = pool.map(populate_runner, range(1, max_k+1))

    hit_runner = functools.partial(count_hits, target_genome=target_genome)
    total_hits_list = pool.map(hit_runner, ((k, bloom) for k, bloom in enumerate(bloom_models, start=1)))

    print(total_hits_list)
    
    return sum(total_hits_list) / len(target_genome)