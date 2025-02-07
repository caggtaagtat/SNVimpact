cd /pathTo/haploSNV/

#For every haplotype csv
while read SAMPLE; do

## Execute R script to introduce sequence variations into the reference sequence
R --vanilla --no-restore --file="/pathTo/haploSNV/Introduce_variations.R" --args $SAMPLE

#End loop
done < /pathTo/haploSNV/xaa

