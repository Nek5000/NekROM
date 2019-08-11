if [[ $1 =~ COPT ]]; then
    ifcopt=1
else
    ifcopt=0
fi

if [[ $1 =~ _L2_ ]]; then
    ifl2='.TRUE.'
    $MOR_DIR/bin/gops cyl_rect_l2
    $MOR_DIR/bin/gbas cyl_rect_l2
else
    ifl2='.FALSE.'
    $MOR_DIR/bin/gops cyl_rect_h10
    $MOR_DIR/bin/gbas cyl_rect_h10
fi

mkdir ops

name="$(echo $1 | perl -pe 's/_(L2|H10)_/_/g')(${ifl2})"
$MOR_DIR/tests/test_template.sh rom_update

ls ../../data/cyl_rect/cyl0.f* > file.list
ls ../../data/cyl_rect/avgcyl0.f* > avg.list
ls bas/bascyl0.f* > bas.list

cp ../../data/cyl_rect/cyl0.f01000 r0.f00001

$SOURCE_ROOT/bin/makenek test

$SOURCE_ROOT/bin/genmap << Z
test
.01
Z

./nek5000 | tee test.log.1

if [[ $ifcopt == 1 ]]; then
   sed -i "s/.*p170/   1.00000     p170/" test.rea 
   ../../tests/tbox
else
   ../../tests/tdragx
fi
