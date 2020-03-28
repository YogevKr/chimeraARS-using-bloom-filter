# chimeraARS using bloom filter

Calculating chimeraARS, using a low space complexity method.

## Installation

```bash
git clone https://github.com/YogevKr/chimeraARS-using-bloom-filter
cd chimeraARS-using-bloom-filter
pip install -r requirements.txt
```

## Usage

```bash
python -m bloomer --host_genome <host_genome path> --target_genome <target_genome path>
```

## OPTIONS
```bash
usage: bloomer.py [-h] --host_genome HOST_GENOME --target_genome TARGET_GENOME
                  [--max_k MAX_K] [--error_rate ERROR_RATE]
                  [--num_processors NUM_PROCESSORS]

optional arguments:
  -h, --help            show this help message and exit
  --host_genome HOST_GENOME
                        Path of host genome FASTA file
  --target_genome TARGET_GENOME
                        Path of target genome FASTA file (Single sequence)
  --max_k MAX_K         Longest sequence
  --error_rate ERROR_RATE
                        Error rate
  --num_processors NUM_PROCESSORS
                        Number of processors
```