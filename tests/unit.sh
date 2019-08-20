if [ "$IPS" = "L2" ]; then
    $MOR_DIR/bin/gops cyl_rect_l2
    $MOR_DIR/bin/gbas cyl_rect_l2
elif [ "$IPS" = "H10" ]; then
    $MOR_DIR/bin/gops cyl_rect_h10
    $MOR_DIR/bin/gbas cyl_rect_h10
else
    echo "inner-product space $IPS not supported..."
    exit(1)
fi

mkdir ops

$MOR_DIR/tests/test_template.sh "$TEST_$IPS_$TYPE"

ls ../../data/cyl_rect/cyl0.f* > file.list
ls ../../data/cyl_rect/avgcyl0.f* > avg.list
ls bas/bascyl0.f* > bas.list

cp ../../data/cyl_rect/cyl0.f01000 r0.f00001

$SOURCE_ROOT/bin/makenek test
$SOURCE_ROOT/bin/genmap << Z
test
.01
Z

./nek5000
