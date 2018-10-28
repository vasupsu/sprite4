#!/bin/bash

PREFIX=$PWD
echo "prefix:" $PREFIX
mkdir -p $PREFIX/bin
builddir=$PREFIX/bin/sprite4_strelka2_install
srcdir=$PREFIX/sprite4_strelka2_modified
mkdir -p $builddir

MPICC=`which mpicc || true`
if test -z $MPICC
then
	echo "No-MPI compilation"
	export CC=gcc
	export CXX=g++
else
	echo "MPI compilation"
	export CC=mpicc
	export CXX=mpicxx
	export MPIFLAG=' -DUSE_MPI'
	export CXXFLAGS=' -DUSE_MPI'
fi

if [ "$(uname)" == "Darwin" ]; then
	echo "Compiling sprite4 for OSX"
	export ECFLAG=' -Xpreprocessor'
	export ELFLAG=' -lomp'
	export CXXFLAGS=' -stdlib=libc++'
fi

#Copy scripts to bin folder
cp sprite4 sprite4-test sprite4-parsnip-test sprite4-generic-test $PREFIX/bin

echo "Compiling minimap2"
cd sprite4_minimap2_modified; make; cp sprite4-minimap2 genFastqIdx $PREFIX/bin; cd ..

echo "Compiling sampa"
$CC -o $PREFIX/bin/sampa ${MPIFLAG} -DUSE_OMP sampa.c ${ECFLAG} -fopenmp ${ELFLAG}

echo "Compiling Strelka2"
cd $builddir
${srcdir}/configure --jobs=4 --prefix=$builddir
make -j4 install
rm -r bootstrap
cd ../..
$CC -o $PREFIX/bin/bamHeaderFile bamHeaderFile.c -I${builddir}/redist/htslib-1.5-1-g03b35d2/htslib -L${builddir}/redist/htslib-1.5-1-g03b35d2  -lhts -lz  -lbz2 -llzma -lpthread -ldeflate
echo "Compiling parsnip"
$CC -o $PREFIX/bin/parsnip parsnip.c -I${builddir}/redist/htslib-1.5-1-g03b35d2/htslib -L${builddir}/redist/htslib-1.5-1-g03b35d2 -lhts -lz  -lbz2 -llzma -lpthread -ldeflate ${MPIFLAG} -DUSE_OMP ${ECFLAG} -fopenmp ${ELFLAG}

echo "Compiling mergevcf"
$CC -o $PREFIX/bin/mergevcf mergeVCF.c


ln -s $builddir/libexec/GetChromDepth $builddir/libexec/starling2 $PREFIX/bin
