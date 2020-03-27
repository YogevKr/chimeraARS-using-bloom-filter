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


def chimera_ars_score(host_genome, target_genome):
    total_hits_list = []

    for k in range(1,50):
        error_rate = 0.01
        max_elements = 4**k if 4**k < len(host_genome) else len(host_genome)

        bloom = BloomFilter(max_elements=max_elements, error_rate=error_rate)
        print(f"k={k}, max_elements={max_elements}, num_bits_m={bloom.num_bits_m}")

        for s in window(host_genome, n=k):
            bloom.add("".join(s))

        k_hits= sum(map(lambda x: x in bloom, window(target_genome, n=k)))

        total_hits_list.append(k_hits)
        print(f"k_hits={k_hits}")

        if k_hits <= 1: break

    print(total_hits_list)
    
    return sum(total_hits_list) / len(target_genome)

