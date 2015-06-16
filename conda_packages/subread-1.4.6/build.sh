
cd src
make -f Makefile.Linux
cd ..

#organize
cp -R bin $PREFIX/

mkdir -p $PREFIX/data/subread
cp -R annotation $PREFIX/data/subread

mkdir -p $PREFIX/share/examples/subread
cp -R test/* $PREFIX/share/examples/subread
cp -R doc/* $PREFIX/share/examples/subread
