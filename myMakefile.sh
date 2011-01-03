#!/bin/bash

CC="ccache g++ -I/opt/intel/ipp/6.1.2.051/ia32/include -Iinclude -Isrc -I/informatik/home/flick/include -I/usr/include -c -fmessage-length=0"
LD="g++ -L/usr/lib -L/opt/intel/ipp/6.1.2.051/ia32/lib -shared"

for i in $(find src -name *.c++)
do
   echo $i
   ${CC} -O0 -g3 -Wall $i -o Debug/`echo $i | sed s/.c++//`.o
done

${LD} $(find Debug -name *.o) -o Debug/libCommonDataAnalysis.so
