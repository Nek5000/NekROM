cd $MOR_DIR/tests
source travis_fn
mkdir t
cd t
$MOR_DIR/bin/linkm

sed 's/(ltr=.*)/(ltr=80)/g' LMOR > LMOR.tmp
mv LMOR.tmp LMOR

cp $MOR_DIR/tests/test.rea .
cp $MOR_DIR/tests/test.usr .

echo 'test'     > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME

fold_start gsnaps "Get Snapshots"; $MOR_DIR/bin/gsnaps cyl_rect; fold_end gsnaps

if [ "$IPS" = "L2" ]; then
    fold_start gops "Get operators"; $MOR_DIR/bin/gops cyl_rect_l2; fold_end gops
    fold_start gbas "Get basis";     $MOR_DIR/bin/gbas cyl_rect_l2; fold_end gbas
   if [ "$TEST" = "CP" ]; then
       fold_start gcp "Get cp operators"; $MOR_DIR/bin/gcp cyl_rect_l2; fold_end gcp
   elif [ "$TEST" = "CP_SKEW" ]; then
       fold_start gcp "Get cp operators"; $MOR_DIR/bin/gcp cyl_rect_l2_skew; fold_end gcp
   fi
elif [ "$IPS" = "H10" ]; then
    fold_start gops "Get operators"; $MOR_DIR/bin/gops cyl_rect_h10; fold_end gops
    fold_start gbas "Get basis";     $MOR_DIR/bin/gbas cyl_rect_h10; fold_end gbas
   if [ "$TEST" = "CP" ]; then
       fold_start gcp "Get cp operators"; $MOR_DIR/bin/gcp cyl_rect_h10; fold_end gcp
   elif [ "$TEST" = "CP_SKEW" ]; then
       fold_start gcp "Get cp operators"; $MOR_DIR/bin/gcp cyl_rect_h10_skew; fold_end gcp
   fi
else
    echo "inner-product space $IPS not supported..."
    exit 1
fi

mv ops tops
mkdir ops &> /dev/null


if [ "$TYPE" = "UNIT" ]; then
    bash $MOR_DIR/tests/unit.sh
elif [ "$TYPE" = "INTEG" ]; then
    bash $MOR_DIR/tests/integ.sh
else
    echo "type $TYPE not supported..."
    exit 1
fi
