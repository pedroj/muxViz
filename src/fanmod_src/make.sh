#!/bin/sh

echo --- NAUTY ---

cd nauty
./configure
make clean
make

#rm *.o
#gcc -c -O3 *.c

cd ..

echo --- COMPILE ---

rm *.o
g++ -c graph64.cpp 
g++ -c output.cpp 
g++ -c random.cpp 
g++ -c maingraph.cpp 
g++ -c main.cpp  

echo --- LINK ---

g++ -o fanmod graph64.o output.o random.o maingraph.o main.o  nauty/gtools.o nauty/nautil.o nauty/nauty.o nauty/naugraph.o nauty/nautinv.o nauty/schreier.o nauty/naurng.o

strip fanmod

#echo --- ZIP ---

#rm fanmod-linux-gtk-i386.zip

#zip fanmod-linux-gtk-i386.zip fanmod

#scp fanmod-linux-gtk-i386.zip login.minet.uni-jena.de:fanmod-linux-gtk-i386.zip
