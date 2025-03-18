outdir=/srv/2024_05_UPF1_PacBio
experiment_file=$outdir/experiment.txt

let p=0
while read -r f1; do
  samples[p]="${f1}"
  ((++p))
done < $experiment_file

for i in "${samples[@]}"; do
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo ""
  echo "$mytime Start minimap2 mapping of sample ${i}"
  echo ""
  
	/home/volker/Tools/minimap2/minimap2 -t 15 -ax splice:hq -C5 --junc-bed /home/volker/reference/Gencode/gencode.v42.SIRVomeERCCome.annotation_forMinimap2.bed /home/volker/reference/Gencode/GRCh38.primary.SIRVomeERCCome.fa ${outdir}/FASTQ/${i}.fastq.gz > ${outdir}/BAM/${i}.sam
	
	mytime=$(date "+%Y-%m-%d %H:%M:%S")
	echo "$mytime Start converting to .bam of sample ${i}"
	  
	samtools view -@ 15 -b ${outdir}/BAM/${i}.sam -o ${outdir}/BAM/${i}.bam
	
	mytime=$(date "+%Y-%m-%d %H:%M:%S")
	echo "$mytime Start sorting sample ${i}"

	samtools sort ${outdir}/BAM/${i}.bam -o ${outdir}/BAM/${i}.sort.bam
	
	mytime=$(date "+%Y-%m-%d %H:%M:%S")
	echo "$mytime Start indexing sample ${i}"

	samtools index ${outdir}/BAM/${i}.sort.bam ${outdir}/BAM/${i}.sort.bam.bai
  
done
