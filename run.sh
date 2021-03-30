# $1 cell numbers along each direction
# $2 Re number
./main $1 $2
mv X.dat Y.dat u.dat v.dat ./verify
cd verify
./plot.sh $2 "Re_$2_$1.png"