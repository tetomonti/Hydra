
make
make razip
make bgzip
mkdir $PREFIX/bin

cp -a samtools $PREFIX/bin/
cp -a razip $PREFIX/bin/
cp -a bgzip $PREFIX/bin/

cp -a bcftools/bcftools $PREFIX/bin/
#cp -a bcftools/bcf-fix.pl $PREFIX/bin/
#cp -a bcftools/vcfutils.pl $PREFIX/bin/

make dylib

mkdir $PREFIX/lib
cp -a libbam.so.1 $PREFIX/lib/

make lib
cp -a libbam.a $PREFIX/lib/

mkdir -p $PREFIX/include/bam/
cp -a *.h $PREFIX/include/bam/

mkdir -p $PREFIX/man/man1/
cp -a samtools.1 $PREFIX/man/man1/

cp -a misc/blast2sam.pl $PREFIX/bin/
cp -a misc/bowtie2sam.pl $PREFIX/bin/
cp -a misc/export2sam.pl $PREFIX/bin/
cp -a misc/interpolate_sam.pl $PREFIX/bin/
cp -a misc/maq2sam-long $PREFIX/bin/
cp -a misc/maq2sam-short $PREFIX/bin/
cp -a misc/md5fa $PREFIX/bin/
cp -a misc/md5sum-lite $PREFIX/bin/
cp -a misc/novo2sam.pl $PREFIX/bin/
cp -a misc/psl2sam.pl $PREFIX/bin/
cp -a misc/sam2vcf.pl $PREFIX/bin/
cp -a misc/samtools.pl $PREFIX/bin/
cp -a misc/soap2sam.pl $PREFIX/bin/
cp -a misc/varfilter.py $PREFIX/bin/
cp -a misc/wgsim $PREFIX/bin/
cp -a misc/wgsim_eval.pl $PREFIX/bin/
cp -a misc/zoom2sam.pl $PREFIX/bin/



