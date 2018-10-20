cp $ROOT_DIR/tests/grammian_test.f t.f
$SOURCE_ROOT/bin/makenek test
$SOURCE_ROOT/bin/genmap << Z
baf
.01
Z
./nek5000
