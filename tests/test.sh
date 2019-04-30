cd $MOR_DIR/tests
mkdir t
cd t
$MOR_DIR/bin/linkm

cp $MOR_DIR/tests/test.rea .
cp $MOR_DIR/tests/test.usr .

echo 'test'     > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME

$MOR_DIR/bin/gsnaps cyl_rect

if [ ${TEST: -4} == "UNIT" ]; then
    bash $MOR_DIR/tests/unit.sh $TEST
else
    bash $MOR_DIR/tests/integ.sh $TEST
fi
