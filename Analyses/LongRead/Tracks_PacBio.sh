mydir=/srv/2024_05_UPF1_PacBio
experiment_file=$mydir/experiment.txt
reference="gencode.v42.SIRVomeERCCome"

let p=0
while read -r f1; do
  samples[p]="${f1}"
  ((++p))
done < $experiment_file

 # Generate Tracks folder for easy visualization
if [ ! -d "${mydir}/Tracks/" ]; then
  mkdir ${mydir}/Tracks/
fi

for i in "${samples[@]}"; do

  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo ""
  echo "$mytime Start track generation of sample ${i}"
  echo ""
  
 # Generate bedGraph using Alfred and convert to tdf using IGVtools
  if [ -e "$mydir/Tracks/${i}.tdf" ]; then
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime ${i} BAM file already converted to track"
  else
    if [ ! -d "${mydir}/Tracks/BedGraph/" ]; then
      mytime=$(date "+%Y-%m-%d %H:%M:%S")
      echo "$mytime Make directory ${mydir}/Tracks/BedGraph/"
      mkdir ${mydir}/Tracks/BedGraph/
  fi

    # Use Alfred to generate bedGraphs
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime ${i} create BedGraph"
    alfred tracks -r 1 -o $mydir/Tracks/BedGraph/${i}.bedGraph.gz $mydir/BAM/${i}.sort.bam
    cd /home/volker/Tools/IGV_2.14.1/

    # Convert to tdf
    mytime=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$mytime ${i} convert to TDF"
    /home/volker/Tools/IGV_2.14.1/igvtools totdf $mydir/Tracks/BedGraph/${i}.bedGraph.gz $mydir/Tracks/${i}.tdf $reference
    # Remove bedGraph files to save disk space
    # rm $mydir/Tracks/BedGraph/${i}.bedGraph.gz
    fi
done
