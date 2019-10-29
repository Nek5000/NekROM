cd $MOR_DIR/tests
source travis_fn
mkdir t
cd t
$MOR_DIR/bin/linkm

cp $MOR_DIR/tests/test.rea .
cp $MOR_DIR/tests/test.usr .

echo 'test'     > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME

$MOR_DIR/bin/gsnaps cyl_rect

if [ "$IPS" = "L2" ]; then
    fold_start gops "Get operators"; $MOR_DIR/bin/gops cyl_rect_l2; fold_end gops
    fold_start gbas "Get basis";     $MOR_DIR/bin/gbas cyl_rect_l2; fold_end gbas
elif [ "$IPS" = "H10" ]; then
    fold_start gops "Get operators"; $MOR_DIR/bin/gops cyl_rect_h10; fold_end gops
    fold_start gbas "Get basis";     $MOR_DIR/bin/gbas cyl_rect_h10; fold_end gbas
else
    echo "inner-product space $IPS not supported..."
    exit 1
fi

mkdir ops

if [ "$TYPE" = "UNIT" ]; then
    $MOR_DIR/tests/unit.sh
elif [ "$TYPE" = "INTEG" ]; then
    $MOR_DIR/tests/integ.sh
else
    echo "type $TYPE not supported..."
    exit 1
fi
