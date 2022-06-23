cd $MOR_DIR/tests
source travis_fn
mkdir t
cd t
$MOR_DIR/bin/linkm

sed 's/(ltr=.*)/(ltr=80)/g' LMOR > LMOR.tmp
mv LMOR.tmp LMOR
sed -i.bu "s/lb=.*)/lb=3)/g" LMOR
sed -i.bu "s/ls=.*)/ls=3)/g" LMOR
sed -i.bu "s/ledvis=.*)/ledvis=1)/g" LMOR
sed -i.bu "s/leleg=.*)/leleg=48)/g" SIZE
sed -i.bu "s/ldimt=.*)/ldimt=4)/g" SIZE

cp $MOR_DIR/tests/rbf_unit.rea .
cp $MOR_DIR/tests/rbf_unit.usr .
cp $MOR_DIR/tests/rbf_unit.mor.

echo 'rbf_unit'     > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME

fold_start grbf "Get rbf weights"; $MOR_DIR/bin/grbf rbf; fold_end grbf

$MOR_DIR/tests/test_template.sh "${TEST}_${TYPE}"

fold_start makenek Makenek
$SOURCE_ROOT/bin/makenek rbf_unit
fold_end makenek

fold_start genmap Genmap
$SOURCE_ROOT/bin/genmap << Z
rbf_unit
.01
Z
fold_end genmap

echo 99 > ecode
 
./nek5000 | tee logfile
iexit=$(cat ./ecode)

if [ "$iexit" != "0" ]; then cp logfile fail.log; fi
exit $iexit
