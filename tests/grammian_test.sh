cp $ROOT_DIR/tests/grammian_test.f t.f
test_name=grammian_test $ROOT_DIR/tests/test_template.sh
ls $ROOT_DIR/data/baf | sed 's/^/..\/..\/data\/baf\//g' > file.list
$SOURCE_ROOT/bin/makenek test
$SOURCE_ROOT/bin/genmap << Z
baf
.01
Z
./nek5000
