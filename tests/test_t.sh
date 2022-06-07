cd $MOR_DIR/tests
source travis_fn
mkdir t
cd t
$MOR_DIR/bin/linkm

sed 's/(ltr=.*)/(ltr=80)/g' LMOR > LMOR.tmp
mv LMOR.tmp LMOR

cp $MOR_DIR/tests/test1.rea .
cp $MOR_DIR/tests/test1.mor .
cp $MOR_DIR/tests/test1.usr .

echo 'test1'     > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME

fold_start gsnaps "Get Snapshots"; $MOR_DIR/bin/gsnaps ann; fold_end gsnaps

mkdir ops &> /dev/null

$MOR_DIR/tests/test_template.sh rom_update

ls ../../data/ann/ann0.f* > file.list
ls ../../data/ann/ann0.f00002 > ic.list
cp ../../data/ann/bas0.f00001 .

sed -i.bu "s/lb=.*)/lb=2)/g" LMOR
sed -i.bu "s/ls=.*)/ls=2)/g" LMOR
sed -i.bu "s/lx1=.*)/lx1=12)/g" SIZE
sed -i.bu "s/lxd=.*)/lxd=18)/g" SIZE

$SOURCE_ROOT/bin/makenek test1
$SOURCE_ROOT/bin/genmap << Z
test1
.01
Z

echo '700000' > Grashof

mpiexec -np 1 ./nek5000 | tee test.log 
cp test.log logfile

echo 99 > ecode

mpiexec -np $NP ./nek5000 | tee test.log | grep -v 'tmax'

cp test.log logfile

grep 'tmax' test.log > nu.log
head nu.log
tail nu.log

../$SCR
iexit=$(cat ./ecode)
if [ "$iexit" != "0" ]; then cp logfile fail.log; fi

ls -latr

exit $iexit
