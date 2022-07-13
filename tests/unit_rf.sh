cd $MOR_DIR/tests
source travis_fn
mkdir t
cd t
$MOR_DIR/bin/linkm

sed 's/(ltr=.*)/(ltr=80)/g' LMOR > LMOR.tmp
mv LMOR.tmp LMOR

cp $MOR_DIR/tests/rf_unit.rea .
cp $MOR_DIR/tests/rf_unit.mor .
cp $MOR_DIR/tests/rf_unit.usr .

echo 'rf_unit'     > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME

$MOR_DIR/bin/gsnaps rft
ls ../../data/rft/chan0.f0000* > file.list

$MOR_DIR/tests/test_template.sh "${TEST}_${TYPE}"

sed -i.bu "s/lb=.*)/lb=2)/g" LMOR
sed -i.bu "s/ls=.*)/ls=2)/g" LMOR
sed -i.bu "s/lx2=.*)/lx2=lx1)/g" SIZE
sed -i.bu "s/leleg=.*)/leleg=48)/g" SIZE
sed -i.bu "s/ldimt=.*)/ldimt=4)/g" SIZE
fold_start makenek Makenek
$SOURCE_ROOT/bin/makenek rf_unit
fold_end makenek

fold_start genmap Genmap
$SOURCE_ROOT/bin/genmap << Z
rf_unit
.01
Z
fold_end genmap

./nek5000 | tee logfile
$MOR_DIR/bin/diff_ascii
iexit=$(cat ./ecode)

if [ "$iexit" != "0" ]; then cp logfile fail.log; fi
exit $iexit
