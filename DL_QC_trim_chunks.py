def write_DL_QC_chunks_paired_separate(SRR_list_file, out_dir, output_bash, group_size = 8, genome = 'mm10_100', SRR_list_server_fileName = 'SRR_list'):

    #this method is for LPC, creates separate download and then processing scripts since can't do both at same time
    fw1 = open(output_bash+"_DL", 'w')
    fw = open(output_bash+"_process",'w')
    
    fw1.write(
'''cd %s
for i in `cat %s`;
do
    echo "wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${i:0: 6}/$i/$i.sra" >>DL_only_$i
done\n''' % (out_dir, SRR_list_server_fileName))
    
    srr_list = open(SRR_list_file, 'r').read().splitlines()
    chunk_list = [srr_list[i:i + group_size] for i in range(0, len(srr_list), group_size)]
    
    for n in range(len(chunk_list)):
        fw1.write('#Processing chunk %s\n'% n)
        arr_string = ''
        chunk = chunk_list[n]
        for samp in chunk:
            fw1.write('bash DL_only_%s &\n'% samp)
            arr_string += '"%s" ' % samp
        fw1.write('wait\n')
    fw.write(
'''cd %s
module add STAR
module add fastqc
module add sratoolkit
module add jdk
module add python/3.6.1
source /home/kyang1/env/bin/activate
mkdir STAR
for i in `cat %s`;
do
    #echo "wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${i:0: 6}/$i/$i.sra" >>DL_and_process_$i
    echo "fastq-dump ./$i.sra --split-3" >>DL_and_process_$i
    echo "rm $i.sra" >> DL_and_process_$i
    echo "mkdir $i\_trimtest" >> DL_and_process_$i
    echo "mv $i\_1.fastq ./$i\_trimtest/" >> DL_and_process_$i
    echo "mv $i\_2.fastq ./$i\_trimtest/" >> DL_and_process_$i
    echo "cd $i\_trimtest" >> DL_and_process_$i
    echo "mkdir ./fastqc_1" >> DL_and_process_$i
    echo "mkdir ./fastqc_2" >> DL_and_process_$i
    echo "fastqc $i\_1.fastq -o ./fastqc_1 &" >> DL_and_process_$i
    echo "fastqc $i\_2.fastq -o ./fastqc_2" >> DL_and_process_$i
    echo "wait" >> DL_and_process_$i
    echo "/home/kyang1/TrimGalore-0.4.5/trim_galore --paired -stringency 5 -length 35 -q 20 -o ./ $i\_1.fastq $i\_2.fastq" >> DL_and_process_$i
    echo "mkdir ./fastqc_1_trimmed" >> DL_and_process_$i
    echo "mkdir ./fastqc_2_trimmed" >> DL_and_process_$i
    echo "fastqc $i\_1_val_1.fq -o ./fastqc_1_trimmed &" >> DL_and_process_$i
    echo "fastqc $i\_2_val_2.fq -o ./fastqc_2_trimmed" >> DL_and_process_$i
    echo "wait" >> DL_and_process_$i
    echo "rm ./*.fastq" >> DL_and_process_$i
    echo "cd .." >> DL_and_process_$i
done\n''' % (out_dir, SRR_list_server_fileName))
    

    
    for n in range(len(chunk_list)):
        fw.write('#Processing chunk %s\n'% n)
        arr_string = ''
        chunk = chunk_list[n]
        for samp in chunk:
            fw.write('bash DL_and_process_%s &\n'% samp)
            arr_string += '"%s" ' % samp
        fw.write('wait\n')
        fw.write('declare -a arr=(%s)\n' % arr_string.strip())
        fw.write(
'''
for i in "${arr[@]}";
do
    mkdir ./STAR/$i
    cd ./STAR/$i
    STAR --genomeDir /project/barash_hdr1/STAR_genomes/%s/ --readFilesIn ../../$i\_trimtest/$i\_1_val_1.fq ../../$i\_trimtest/$i\_2_val_2.fq --runThreadN 7 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./$i. --outSAMattributes All --alignSJoverhangMin 8
    rm ../../$i\_trimtest/*.fq
    cd ../../
done\n''' % genome)
    fw.close()
    fw1.close()

def write_DL_QC_chunks_single_separate(SRR_list_file, out_dir, output_bash, group_size = 8, genome = 'dr11_100', SRR_list_server_fileName = 'SRR_list'):

    #this method is for LPC, creates separate download and then processing scripts since can't do both at same time
    fw1 = open(output_bash+"_DL", 'w')
    fw = open(output_bash+"_process",'w')
    
    fw1.write(
'''cd %s
for i in `cat %s`;
do
    echo "wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${i:0: 6}/$i/$i.sra" >>DL_only_$i
done\n''' % (out_dir, SRR_list_server_fileName))

    srr_list = open(SRR_list_file, 'r').read().splitlines()
    chunk_list = [srr_list[i:i + group_size] for i in range(0, len(srr_list), group_size)]

    for n in range(len(chunk_list)):
        fw1.write('#Processing chunk %s\n'% n)
        arr_string = ''
        chunk = chunk_list[n]
        for samp in chunk:
            fw1.write('bash DL_only_%s &\n'% samp)
            arr_string += '"%s" ' % samp
        fw1.write('wait\n')
    fw.write(
'''cd %s
module add STAR
module add fastqc
module add sratoolkit
module add jdk
module add python/3.6.1
source /home/kyang1/env/bin/activate
mkdir STAR
for i in `cat %s`;
do
    #echo "wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${i:0: 6}/$i/$i.sra" >>DL_and_process_$i
    echo "fastq-dump ./$i.sra --split-3" >>DL_and_process_$i
    echo "rm $i.sra" >> DL_and_process_$i
    echo "mkdir $i\_trimtest" >> DL_and_process_$i
    echo "mv $i.fastq ./$i\_trimtest/" >> DL_and_process_$i
    #echo "mv $i\_2.fastq ./$i\_trimtest/" >> DL_and_process_$i
    echo "cd $i\_trimtest" >> DL_and_process_$i
    echo "mkdir ./fastqc_1" >> DL_and_process_$i
    #echo "mkdir ./fastqc_2" >> DL_and_process_$i
    echo "fastqc $i.fastq -o ./fastqc_1 &" >> DL_and_process_$i
    #echo "fastqc $i\_2.fastq -o ./fastqc_2" >> DL_and_process_$i
    echo "wait" >> DL_and_process_$i
    echo "/home/kyang1/TrimGalore-0.4.5/trim_galore -stringency 5 -length 35 -q 20 -o ./ $i.fastq" >> DL_and_process_$i
    echo "mkdir ./fastqc_1_trimmed" >> DL_and_process_$i
    #echo "mkdir ./fastqc_2_trimmed" >> DL_and_process_$i
    echo "fastqc $i\_trimmed.fq -o ./fastqc_1_trimmed &" >> DL_and_process_$i
    #echo "fastqc $i\_2_val_2.fq -o ./fastqc_2_trimmed" >> DL_and_process_$i
    echo "wait" >> DL_and_process_$i
    echo "rm ./*.fastq" >> DL_and_process_$i
    echo "cd .." >> DL_and_process_$i
done\n''' % (out_dir, SRR_list_server_fileName))
    
    
    for n in range(len(chunk_list)):
        fw.write('#Processing chunk %s\n'% n)
        arr_string = ''
        chunk = chunk_list[n]
        for samp in chunk:
            fw.write('bash DL_and_process_%s &\n'% samp)
            arr_string += '"%s" ' % samp
        fw.write('wait\n')
        fw.write('declare -a arr=(%s)\n' % arr_string.strip())
        fw.write(
'''
for i in "${arr[@]}";
do
    mkdir ./STAR/$i
    cd ./STAR/$i
    STAR --genomeDir /project/barash_hdr1/STAR_genomes/%s/ --readFilesIn ../../$i\_trimtest/$i\_trimmed.fq --runThreadN 7 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./$i. --outSAMattributes All --alignSJoverhangMin 8
    #rm ../../$i\_trimtest/*.fq
    cd ../../
done\n''' % genome)
    fw.close()
    fw1.close()


def write_DL_QC_trim_chunks_single(SRR_list_file, out_dir, output_bash, group_size = 8, genome = 'mm10', SRR_list_server_fileName = 'SRR_list'):

    fw = open(output_bash, 'w')
    
    fw.write(
'''cd %s
mkdir STAR
for i in `cat %s`;
do
    echo "wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${i:0: 6}/$i/$i.sra" >>DL_and_process_$i
    echo "fastq-dump ./$i.sra --split-3" >>DL_and_process_$i
    echo "rm $i.sra" >> DL_and_process_$i
    echo "mkdir $i\_trimtest" >> DL_and_process_$i
    echo "mv $i.fastq ./$i\_trimtest/" >> DL_and_process_$i
    #echo "mv $i\_2.fastq ./$i\_trimtest/" >> DL_and_process_$i
    echo "cd $i\_trimtest" >> DL_and_process_$i
    echo "mkdir ./fastqc_1" >> DL_and_process_$i
    #echo "mkdir ./fastqc_2" >> DL_and_process_$i
    echo "fastqc $i.fastq -o ./fastqc_1 &" >> DL_and_process_$i
    #echo "fastqc $i\_2.fastq -o ./fastqc_2" >> DL_and_process_$i
    echo "wait" >> DL_and_process_$i
    echo "trim_galore -stringency 5 -length 35 -q 20 -o ./ $i.fastq" >> DL_and_process_$i
    echo "mkdir ./fastqc_1_trimmed" >> DL_and_process_$i
    #echo "mkdir ./fastqc_2_trimmed" >> DL_and_process_$i
    echo "fastqc $i\_trimmed.fq -o ./fastqc_1_trimmed &" >> DL_and_process_$i
    #echo "fastqc $i\_2_val_2.fq -o ./fastqc_2_trimmed" >> DL_and_process_$i
    echo "wait" >> DL_and_process_$i
    echo "rm ./*.fastq" >> DL_and_process_$i
    echo "cd .." >> DL_and_process_$i
done\n''' % (out_dir, SRR_list_server_fileName))
    
    srr_list = open(SRR_list_file, 'r').read().splitlines()
    chunk_list = [srr_list[i:i + group_size] for i in range(0, len(srr_list), group_size)]
    
    for n in range(len(chunk_list)):
        fw.write('#Processing chunk %s\n'% n)
        arr_string = ''
        chunk = chunk_list[n]
        for samp in chunk:
            fw.write('bash DL_and_process_%s &\n'% samp)
            arr_string += '"%s" ' % samp
        fw.write('wait\n')
        fw.write('declare -a arr=(%s)\n' % arr_string.strip())
        fw.write(
'''
for i in "${arr[@]}";
do
    mkdir ./STAR/$i
    cd ./STAR/$i
    STAR --genomeDir /data/STAR_DATA/%s/ --readFilesIn ../../$i\_trimtest/$i\_trimmed.fq --runThreadN 16 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./$i. --outSAMattributes All --alignSJoverhangMin 8
    #rm ../../$i\_trimtest/*.fq
    cd ../../
done\n''' % genome)
    fw.close()


def write_DL_QC_trim_chunks(SRR_list_file, out_dir, output_bash, group_size = 8, genome = 'mm10', SRR_list_server_fileName = 'SRR_list'):

    fw = open(output_bash, 'w')
    
    fw.write(
'''cd %s
mkdir STAR
for i in `cat %s`;
do
    echo "wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${i:0: 6}/$i/$i.sra" >>DL_and_process_$i
    echo "fastq-dump ./$i.sra --split-3" >>DL_and_process_$i
    echo "rm $i.sra" >> DL_and_process_$i
    echo "mkdir $i\_trimtest" >> DL_and_process_$i
    echo "mv $i\_1.fastq ./$i\_trimtest/" >> DL_and_process_$i
    echo "mv $i\_2.fastq ./$i\_trimtest/" >> DL_and_process_$i
    echo "cd $i\_trimtest" >> DL_and_process_$i
    echo "mkdir ./fastqc_1" >> DL_and_process_$i
    echo "mkdir ./fastqc_2" >> DL_and_process_$i
    echo "fastqc $i\_1.fastq -o ./fastqc_1 &" >> DL_and_process_$i
    echo "fastqc $i\_2.fastq -o ./fastqc_2" >> DL_and_process_$i
    echo "wait" >> DL_and_process_$i
    echo "trim_galore --paired -stringency 5 -length 35 -q 20 -o ./ $i\_1.fastq $i\_2.fastq" >> DL_and_process_$i
    echo "mkdir ./fastqc_1_trimmed" >> DL_and_process_$i
    echo "mkdir ./fastqc_2_trimmed" >> DL_and_process_$i
    echo "fastqc $i\_1_val_1.fq -o ./fastqc_1_trimmed &" >> DL_and_process_$i
    echo "fastqc $i\_2_val_2.fq -o ./fastqc_2_trimmed" >> DL_and_process_$i
    echo "wait" >> DL_and_process_$i
    echo "rm ./*.fastq" >> DL_and_process_$i
    echo "cd .." >> DL_and_process_$i
done\n''' % (out_dir, SRR_list_server_fileName))
    
    srr_list = open(SRR_list_file, 'r').read().splitlines()
    chunk_list = [srr_list[i:i + group_size] for i in range(0, len(srr_list), group_size)]
    
    for n in range(len(chunk_list)):
        fw.write('#Processing chunk %s\n'% n)
        arr_string = ''
        chunk = chunk_list[n]
        for samp in chunk:
            fw.write('bash DL_and_process_%s &\n'% samp)
            arr_string += '"%s" ' % samp
        fw.write('wait\n')
        fw.write('declare -a arr=(%s)\n' % arr_string.strip())
        fw.write(
'''
for i in "${arr[@]}";
do
    mkdir ./STAR/$i
    cd ./STAR/$i
    STAR --genomeDir /data/STAR_DATA/%s/ --readFilesIn ../../$i\_trimtest/$i\_1_val_1.fq ../../$i\_trimtest/$i\_2_val_2.fq --runThreadN 16 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./$i. --outSAMattributes All --alignSJoverhangMin 8
    rm ../../$i\_trimtest/*.fq
    cd ../../
done\n''' % genome)
    fw.close()


def write_unzip_QC_trim_chunks(SRR_list_file, out_dir, output_bash, group_size = 8, genome = 'mm10', SRR_list_server_fileName = 'SRR_list'):

    fw = open(output_bash, 'w')
    
    fw.write(
'''cd %s
mkdir STAR
for i in `cat %s`;
do
    echo "gunzip ./fastq/$i\_R1_001.fastq.gz" >> DL_and_process_$i
    echo "gunzip ./fastq/$i\_R2_001.fastq.gz" >> DL_and_process_$i
    echo "mkdir ./$i\_trimtest/" >> DL_and_process_$i
    echo "mv ./fastq/$i\_R1_001.fastq ./$i\_trimtest/$i\_1.fastq" >> DL_and_process_$i
    echo "mv ./fastq/$i\_R2_001.fastq ./$i\_trimtest/$i\_2.fastq" >> DL_and_process_$i
    echo "cd $i\_trimtest" >> DL_and_process_$i
    echo "mkdir ./fastqc_1" >> DL_and_process_$i
    echo "mkdir ./fastqc_2" >> DL_and_process_$i
    echo "fastqc $i\_1.fastq -o ./fastqc_1 &" >> DL_and_process_$i
    echo "fastqc $i\_2.fastq -o ./fastqc_2" >> DL_and_process_$i
    echo "wait" >> DL_and_process_$i
    echo "trim_galore --paired -stringency 5 -length 35 -q 20 -o ./ $i\_1.fastq $i\_2.fastq" >> DL_and_process_$i
    echo "mkdir ./fastqc_1_trimmed" >> DL_and_process_$i
    echo "mkdir ./fastqc_2_trimmed" >> DL_and_process_$i
    echo "fastqc $i\_1_val_1.fq -o ./fastqc_1_trimmed &" >> DL_and_process_$i
    echo "fastqc $i\_2_val_2.fq -o ./fastqc_2_trimmed" >> DL_and_process_$i
    echo "wait" >> DL_and_process_$i
    echo "rm ./*.fastq" >> DL_and_process_$i
    echo "cd .." >> DL_and_process_$i
done\n''' % (out_dir, SRR_list_server_fileName))
    
    srr_list = open(SRR_list_file, 'r').read().splitlines()
    chunk_list = [srr_list[i:i + group_size] for i in range(0, len(srr_list), group_size)]
    
    for n in range(len(chunk_list)):
        fw.write('#Processing chunk %s\n'% n)
        arr_string = ''
        chunk = chunk_list[n]
        for samp in chunk:
            fw.write('bash DL_and_process_%s &\n'% samp)
            arr_string += '"%s" ' % samp
        fw.write('wait\n')
        fw.write('declare -a arr=(%s)\n' % arr_string.strip())
        fw.write(
'''
for i in "${arr[@]}";
do
    mkdir ./STAR/$i
    cd ./STAR/$i
    STAR --genomeDir /data/STAR_DATA/%s/ --readFilesIn ../../$i\_trimtest/$i\_1_val_1.fq ../../$i\_trimtest/$i\_2_val_2.fq --runThreadN 16 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./$i. --outSAMattributes All --alignSJoverhangMin 8
    rm ../../$i\_trimtest/*.fq
    cd ../../
done\n''' % genome)
    fw.close()
