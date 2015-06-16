#!/bin/bash

./configure \
    --enable-shared \
    --enable-ipv6 \
    --enable-unicode=ucs4 \
    --prefix=$PREFIX
make
make install
  
