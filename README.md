# kenezam

## trimming (GSE137757/fastq_raw_annotated/extra_trim.sh)

~/bin/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 20 ${d} ../fastq_trimmed/${d:0:-11}4s.fastq.gz ILLUMINACLIP:/home/danilb/bin/adapters.fa:2:30:5 SLIDINGWINDOW:4:20 MINLEN:20

## Length distribution (pictures_len_dist.sh)
unzip -p GSE135860/fastq_wo_ribo/GSE135860_ribo_virus_HBV_72h_1_5s_fastqc.zip GSE135860_ribo_virus_HBV_72h_1_5s_fastqc/Images/sequence_length_distribution.png > res2.png

 ~/bin/magick res0.png res1.png res2.png +append ../len_dist_comp/ribo/${e}.png;
 
 ## Picture
 ### мдва
 ![alt text](https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true)
