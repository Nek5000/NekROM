cd $MOR_DIR/tests
source travis_fn
mkdir t
cd t
$MOR_DIR/bin/linkm

cp $MOR_DIR/tests/test.rea .
cp $MOR_DIR/tests/test.usr .

echo 'test'     > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME

fold_start gsnaps "Get Snapshots"; $MOR_DIR/bin/gsnaps cyl_rect; fold_end gsnaps

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

mv ops tops
mkdir ops &> /dev/null

if [ "$TEST" = "CP" ]; then
    fold_start gbas "Get cp operators"; $MOR_DIR/bin/gcp cyl_rect_h10; fold_end gbas
fi

if [ "$TYPE" = "UNIT" ]; then
    bash $MOR_DIR/tests/unit.sh
elif [ "$TYPE" = "INTEG" ]; then
    bash $MOR_DIR/tests/integ.sh
else
    echo "type $TYPE not supported..."
    exit 1
fi
