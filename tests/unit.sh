$MOR_DIR/tests/test_template.sh "${TEST}_${IPS}_${TYPE}"

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
