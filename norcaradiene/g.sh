#!/bin/bash
x=$1
y=${x%.*}

babel $x $y-freq.gjf
z=$y-freq
echo "%nproc=4" > $z.com
echo "# freq mp2/aug-cc-pvdz" >> $z.com
tail -n+2 $z.gjf >> $z.com
rm $z.gjf

/home/eliot/software/g09E01/g09 < $z.com > $z.log &



