# Activate the specified Conda environment
#conda activate VIP2B.1.1

mkdir -p data/

# Download data
if [[ ! -f data/test_seq.R1.fq.gz ]]; then
    wget --user-agent="Mozilla/5.0" -t 3 -O data/test_seq.R1.fq.gz https://figshare.com/ndownloader/files/52717946
fi
if [[ ! -f data/test_seq.R2.fq.gz ]]; then
    wget --user-agent="Mozilla/5.0" -t 3 -O data/test_seq.R2.fq.gz https://figshare.com/ndownloader/files/52717949
fi

# Run pipeline
time python3 ../bin/VIP2B.py -i data.list
