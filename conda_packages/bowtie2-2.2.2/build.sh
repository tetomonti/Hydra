#./configure --prefix $PREFIX
make
#make install
chmod a+x bowtie2
cp bowtie2 \
   bowtie2-align-s \
   bowtie2-align-l \
   bowtie2-build \
   bowtie2-build-s \
   bowtie2-build-l \
   bowtie2-inspect \
   bowtie2-inspect-s \
   bowtie2-inspect-l \
   $PREFIX/bin

chmod a+x scripts/*
cp scripts/* $PREFIX/bin
