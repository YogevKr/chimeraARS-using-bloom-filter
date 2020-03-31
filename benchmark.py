from itertools import islice
from Bio import SeqIO, Seq
from collections import defaultdict
import json

from bloomer import chimera_ars_score, naive_chimera_ars

def load_host_genome(path):
    return sum(
    [g.seq for g in SeqIO.parse(path, "fasta")], Seq.Seq("")
).upper()

def load_target_genome(path):
    return tuple(
    islice(SeqIO.parse(path, "fasta"), 1)
)[0].seq.upper()

host_Ecoli: Seq.Seq = load_host_genome("./GCF_000005845.2_ASM584v2_genomic.fna")
host_Scer: Seq.Seq = load_host_genome("./GCF_000146045.2_R64_genomic.fna")
target: Seq.Seq = load_target_genome("./target.fna")
    
results = defaultdict(lambda: defaultdict())

results["target"]["length"] = len(target)
results["Ecoli"]["length"] = len(host_Ecoli)
results["Scer"]["length"] = len(host_Scer)

results["Ecoli"]["naive_hits"] = naive_chimera_ars(host_Ecoli, target)
results["Scer"]["naive_hits"] = naive_chimera_ars(host_Scer, target)

results["Ecoli"]["estimation_details"] = []
results["Scer"]["estimation_details"] = []

for e in [0.5, 0.4, 0.3, 0.2, 0.1, 0.01, 0.001, 0.0001, 0.00001]:
    print(f"## Error {e} ##")
    ecoli_score, ecoli_hits = chimera_ars_score(
        host_genome=host_Ecoli,
        target_genome=target,
        max_k=20,
        error_rate=e,
        num_processors=0,
    )
        
    results["Ecoli"]["estimation_details"].append((e, ecoli_score, ecoli_hits))
    
    scer_score, scer_hits = chimera_ars_score(
        host_genome=host_Scer,
        target_genome=target,
        max_k=20,
        error_rate=e,
        num_processors=0,
    )
        
    results["Scer"]["estimation_details"].append((e, scer_score, scer_hits))
        
print(dict(results))

with open('results.json', 'w') as outfile:
    json.dump(dict(results), outfile)