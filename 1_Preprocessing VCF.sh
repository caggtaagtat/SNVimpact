## Step 1: Splitting large VCF into VCF per individual
file=hprc-v1.1-mc-chm13.GRCh38.vcfbub.a100k.wave.vcf.gz

for sample in `bcftools query -l $file`; do
  bcftools view -c1 -Oz -s $sample -o ${file/.vcf*/.$sample.small.vcf.gz} $file
done


## Step 2: Convert Individual VCF to hap Format for later Separation of haplotype 1 and 2
find /pathTo/haploSNV/newVCF/ -name "*small*" | while read SMALLFILE

do
#echo $SMALLFILE
f="$(basename -- $SMALLFILE)"
f2="$(basename -- $SMALLFILE ".small.vcf.gz")"
bcftools convert --hapsample --vcf-ids $f -o $f2".hap.vcf"

done


## Step 3: For every hap file, generate SNV csv file per haplotype using custom R script "create_csv.R"

cd /pathTo/haploSNV/newVCF/haps/csv/

#For every file in the given path, showing a "hap" at the end of the filename
find /pathTo/haploSNV/newVCF/haps/ -name "*.hap" | while read SAMPLE

do

## Execute R script to seperate haplotypes
R --vanilla --no-restore --file="/pathTo/haploSNV/2_Create_csv.R" --args $SAMPLE

done




