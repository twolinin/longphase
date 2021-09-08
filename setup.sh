cd htslib
autoreconf -i
./configure
make -j 4
cd ..
autoreconf -i
./configure
make -j 4

