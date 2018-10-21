# MOR Tests

test_unit() {
$ROOT_DIR/tests/test_template.sh $1
ls $ROOT_DIR/data/baf | sed 's/^/..\/..\/data\/baf\//g' > file.list
$SOURCE_ROOT/bin/makenek test
$SOURCE_ROOT/bin/genmap << Z
baf
.01
Z
./nek5000
}

cd $ROOT_DIR/cases/baf
$ROOT_DIR/bin/linkc

echo 'baf'     > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME

$ROOT_DIR/bin/gsnaps baf
$ROOT_DIR/bin/gops   baf

case "$TEST" in
    GRAMMIAN_UNIT)
        test_unit grammian_test
        ;;
    EIGENVECTOR_UNIT)
        test_unit eigenvector_test
        ;;
    BASES_UNIT)
        test_unit bases_test
        ;;
    BAFF_COMP)
        $ROOT_DIR/tests/comp-baff-test
        ;;
    BAF_COMP)
        $ROOT_DIR/tests/comp-baf-test
        ;;
    BAF_INTEG)
        $ROOT_DIR/tests/baf-test
        ;;
    *)
        echo 'did not specify test...'
        ;;
esac
