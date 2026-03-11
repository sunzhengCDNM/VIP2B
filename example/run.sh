# Activate the specified Conda environment
#conda activate VIP2B.1.1

mkdir -p data/

# Download data
if [[ ! -f data/test_seq.R1.fq.gz ]]; then
    wget -t 3 -O data/test_seq.R1.fq.gz https://zenodo.org/records/18947711/files/test_seq.R1.fq
fi
if [[ ! -f data/test_seq.R2.fq.gz ]]; then
    wget -t 3 -O data/test_seq.R2.fq.gz https://zenodo.org/records/18947711/files/test_seq.R2.fq
fi

# Run pipeline
time python3 ../bin/VIP2B.py -i data.list
