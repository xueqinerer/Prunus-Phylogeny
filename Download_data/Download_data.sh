#1.Download the data from NCBI SRA
#install prefetch
#prefetch the current directory(select the number of SRA)
prefetch SRR15882990 ./
#pretch batch download
prefetch –option-file SRR_Acc_List.txt
#prefetch –max-size 40G SRR1615263，通过使用–max-size参数来设置下载的最大大小，这样可以确保您下载的数据在限制范围内
#SRR_Acc_List.txt
## SRR13783530
## SRR23496001
## SRR22430821
## SRR19658861
## SRR18715737
## SRR15112554

#2.Replace the filename from the NCBI SRA database,Obtaining the name of the target species for later analysis

a.single filename
rename 's/SRR12927737.clean/Prunus_pumila/' *.fq or mv SRR12927737.clean.fq Prunus_pumila.fq

b.Multiple filenames(modify_name.sh)
#define the original and new file names readme.txt
filename_list="./readme.txt"
#loop a list include riginal and new file names
while read -r line;
do
#old name means "SRR12927737"
#new name means "redefine new names"
old_name=$(echo "$line" | awk '{print $1}')
new_name=$(echo "$line" | awk '{print $2}')

#check the old name exist
if [ -e "$old_name" ]
then
#rename the filenames
mv "$old_name" "$new_name" 
echo "Renamed $old_name to $new_name"
else
  echo "File $old_name not found"
fi

done < "$filename_list"
##            V1                    V2
## 1  ERR7250757   Amygdalus_amygdalus
## 2 SRR17221713   Amygdalus_davidiana
## 3 SRR15403606  Amygdalus_kansuensis
## 4 SRR15403597        Amygdalus_mira
## 5  SRR1652475   Amygdalus_mongolica
## 6 SRR13261910 Amygdalus_pedunculata


#3.The format change: fastq-dump
a.sra——fastaq
for i in *.sra;
do
echo $i
fastq-dump --split-3 $i
done

b.fasta——fastq
#single data
seqtk seq Pygeum_topengii.fasta -F  "J" > Pygeum_topengii.fq

#mulitiple data
for sample in *.fasta;do
sample_name=$(basename ${sample} .fasta)
echo $sample_name
seqtk seq $sample -F  "J" > $sample.fq
done

#4.Trimmomatic the data: trimmomatic Bolger et al (2014),Clean the data,in order to generate clean data

(1)single data
a.Paired
java -jar /data/xueqin/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 4 -phred33 Prunus_amygdalus_1.fastq SRR15412870_2.fastq Prunus_amygdalus_1.clean.R1.fq.gz Prunus_amygdalus.clean.unpaired.R1.fq.gz Prunus_amygdalus_2.clean.R2.fq.gz Prunus_amygdalus.clean.unpaired.R2.fq.gz ILLUMINACLIP:/data/xueqin/software/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

b.Single
java -jar /home/xueqin/data/software/Trimmomatic-0.38/trimmomatic-0.38.jar SE -threads 8 -phred33 Amygdalus_pedunculata.fastq Amygdalus_pedunculata.fq.gz ILLUMINACLIP:/home/xueqin/data/software/Trimmomatic-0.38/adapters/TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

(2)multiple data
a.paired_end(trimmomatic_paired.sh)
#!/bin/bash
input_dir=/home/xueqin/data/Prunus_data_analyses/Transcriptomic_all/SRA_Transcriptomic/test
output_dir=/home/xueqin/data/Prunus_data_analyses/Transcriptomic_all/SRA_Transcriptomic/test
#specify the parameter and pathway of the Trimmomatic
trimmomatic_path=/home/xueqin/data/software/Trimmomatic-0.38/trimmomatic-0.38.jar
trimmomatic_params="ILLUMINACLIP:/home/xueqin/data/software/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
#deal with the sample_1.fastq
for sample in ${input_dir}/*_1.fastq; do
#extract the sample name(it means the delete the "_1.fastq")
  sample_name=$(basename ${sample} _1.fastq)
#the left is the sample_1.fastq
  left_reads=${sample}
#the right equal sample_2.fastq
  right_reads=${input_dir}/${sample_name}_2.fastq
  output_prefix=${output_dir}/${sample_name}
#Trimmomatic the data
  java -jar ${trimmomatic_path} PE -threads 8 -phred33 ${left_reads} ${right_reads} ${output_prefix}_1.fq.gz ${output_prefix}_1_unpaired.fastq.gz ${output_prefix}_2.fq.gz ${output_prefix}_2_unpaired.fastq.gz ${trimmomatic_params}
done

b.single_end(trimmomatic_single.sh)
#!/bin/bash
input_dir=/home/xueqin/data/Prunus_dada_analyses/Transcriptomic_all/SRA_Transcriptomic/Fastq_Transcritomic/Singer
output_dir=/home/xueqin/data/Prunus_dada_analyses/Transcriptomic_all/SRA_Transcriptomic
trimmomatic_path=/home/xueqin/data/software/Trimmomatic-0.38/trimmomatic-0.38.jar
trimmomatic_params="ILLUMINACLIP:/home/xueqin/data/software/Trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
for sample in ${input_dir}/*.fastq; do
  sample_name=$(basename ${sample} .fastq)
  reads=${sample}
  output_prefix=${output_dir}/${sample_name}
  java -jar ${trimmomatic_path} SE -threads 8 -phred33 ${reads} ${output_prefix}.fq.gz ${trimmomatic_params}
done

#5. Quality control:fastqc
#look the .gz file
fastqc *.fastq.gz
#summarizes quality control results
multiqc ./

#6.Gunzip the clean data,in order to analysis
for i in *gz;
do
echo $i
gunzip $i
done

