gfortran -c ../src-lib/csor.f
gfortran csor.o ../src-driver/test.f -o test
./test
