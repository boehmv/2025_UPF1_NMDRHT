outdir=/srv/2024_06_UPF1_Nanopore
experiment_file=$outdir/experiment.txt

let p=0
while read -r f1 f2; do
  samples[p]="${f1}"
  ((++p))
done < $experiment_file

source /home/volker/miniconda3/bin/activate activate nanoplot 

for i in "${samples[@]}"; do
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo ""
  echo "$mytime Start BAM QC of sample ${i}"
  echo ""
  
  mkdir -p ${outdir}/QC/BAM/${i}
  
  NanoPlot -t 12 --bam ${outdir}/BAM/${i}.sort.bam --huge --info_in_report  -o ${outdir}/QC/BAM/${i}
  
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo ""
  echo "$mytime Start FASTQ QC of sample ${i}"
  echo ""
  
   mkdir -p ${outdir}/QC/FASTQ/${i}

  NanoPlot -t 12 --fastq ${outdir}/FASTQ/${i}.fastq.gz --huge --info_in_report  -o ${outdir}/QC/FASTQ/${i}
  
done

conda deactivate
