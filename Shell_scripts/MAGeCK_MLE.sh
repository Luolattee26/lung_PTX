cd /data/jxwang_data/yjz/Sample/lung_PTX_count


# MAGeCK COUNT
# TEST123
mageck count -l ../Human_GeCKOv2_Library_A_3_mageck.csv -n sgRNA_count_test123 --fastq ../WGC085985-A02_R1.fastq \
../WGC085985-A03_R1.fastq ../WGC085985-A04_R1.fastq \
../WGC085985-A05_R1.fastq ../WGC085985-A06_R1.fastq \
../WGC085985-A07_R1.fastq ../WGC085985-A08_R1.fastq \
../WGC085985-A09_R1.fastq ../WGC085985-A10_R1.fastq \
../WGC085985-A11_R1.fastq ../WGC085985-A12_R1.fastq \
--sample-label test1.D0,test1.D7,test1.T14,test1.D14,test2.D0,test2.D7,test2.T14,test3.D0,test3.D7,test3.T14,test3.D14


# MAGeCK batch effect removal
# R package


# MAGeCK MLE

# mageck mle --count-table raw_count_batch_removed_13.txt \
# --design-matrix designMTX.txt --norm-method control \
# --control-sgrna nonessential_sgrna_list.txt \
# --output-prefix test13.mle

# if use batch removed data,some negetive read counts will appear

cd /data/jxwang_data/yjz/Sample/lung_PTX_mle

# all 123
mageck mle --count-table ../lung_PTX_count/sgRNA_count_test123.count.txt \
--design-matrix designMTX_all_tests.txt --norm-method control \
--control-sgrna nonessential_sgrna_list.txt \
--output-prefix all_test123.mle --permutation-round 10 \
--remove-outliers --threads 4