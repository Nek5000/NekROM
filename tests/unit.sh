$MOR_DIR/tests/test_template.sh "${TEST}_${IPS}_${TYPE}"

ls ../../data/cyl_rect/cyl0.f* > file.list
ls ../../data/cyl_rect/avgcyl0.f* > avg.list
ls bas/bascyl0.f* > bas.list

cp ../../data/cyl_rect/cyl0.f01000 r0.f00001

fold_start makenek Makenek
$SOURCE_ROOT/bin/makenek test
fold_end makenek

fold_start genmap Genmap
$SOURCE_ROOT/bin/genmap << Z
test
.01
Z
fold_end genmap

echo 99 > ecode
 
./nek5000 | tee logfile
iexit=$(cat ./ecode)

if [ "$iexit" != "0" ]; then cp logfile fail.log; fi
exit $iexit
