#!/bin/bash

# GSLINCLUDE=`pkg-config gsl --cflags`
CCDEBUG="ccache g++ ${GSLINCLUDE} -I. -Iinclude -Isrc -I/usr/include -c -fmessage-length=0 -DVERBOSE -DEXTRA_VERBOSE -DDETAIL_VERBOSE"
CC=${CCDEBUG}
LD="g++ -L/usr/lib -shared"

# Verzeichnisbaum für .o-Dateien erzeugen (weiß der Geier, habe ich von eclipse CDT übernommen)
for i in $(find src -type d)
do
   (mkdir Debug/$i 2>&1 >/dev/null) || true
done

# Alle .o-Dateien erzeugen
for i in $(find src -name *.c++)
do
   echo $i
   ${CC} -O0 -g3 -Wall $i -o Debug/`echo $i | sed s/.c++//`.o 
done

# Linken
echo "Linking ..."
${LD} $(find Debug -name *.o) -o Debug/libCommonDataAnalysis.so

