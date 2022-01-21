#!/usr/bin/bash
for i in x*;
do
echo "#!/usr/bin/bash" >Run.$i.sh
echo "mkdir $i.dir;" >>Run.$i.sh
echo "mv $i $i.dir; " >>Run.$i.sh
echo "cd $i.dir" >>Run.$i.sh
echo "\`perl LTR_detector/Script/Solo_LTR_TSD_checker.pl Genome_seq/Temp/1st.round.ltr.X.masked.genome.fasta  $i & \`" >>Run.$i.sh
echo "cd ..;" >>Run.$i.sh

done
