./configure \
      --prefix="$PREFIX" \
      --enable-R-profiling \
      --enable-memory-profiling \
      --enable-R-shlib \
      --enable-BLAS-shlib \
      --enable-byte-compiled-packages \
      --with-readline \
      --with-tcltk \
      --with-cairo \
      --with-libpng \
      --with-jpeglib \
      --with-recommended-packages \
      --with-x

make
make pdf
make info
pushd ./src/nmath/standalone

  make
  make install
  
popd
make install
make install-pdf
make install-info
