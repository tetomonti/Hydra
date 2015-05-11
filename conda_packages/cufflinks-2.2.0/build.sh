#!/bin/bash

# ./configure --with-boost=$PREFIX \
# 	    --with-boost-thread=$PREFIX/lib/libboost_thread.a \
# 	    --with-boost-system=$PREFIX/lib/libboost_system.a \
# 	    --with-boost-serialization=$PREFIX/lib/libboost_serialization.a \
# 	    --with-bam=$PREFIX \
# 	    --with-eigen=$PREFIX/include


./configure --prefix=$PREFIX \
	    --with-boost=$PREFIX \
	    --with-boost-thread=$PREFIX/lib/libboost_thread.a \
	    --with-boost-system=$PREFIX/lib/libboost_system.a \
	    --with-boost-serialization=$PREFIX/lib/libboost_serialization.a \
	    --with-bam=$PREFIX \
	    --with-eigen=$PREFIX/include


#./configure --prefix $PREFIX
make
make install 
# cp src/cufflinks $PREFIX/bin
# cp src/cuffcompare $PREFIX/bin
# cp src/cuffdiff $PREFIX/bin
# cp src/cuffmerge $PREFIX/bin/cuffmerge
# cp src/gffread $PREFIX/bin
# cp src/gtf_to_sam $PREFIX/bin
# cp src/cuffnorm $PREFIX/bin
# cp src/cuffquant $PREFIX/bin
