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
    INITIAL_CONDITION_UNIT)
        test_unit initial_condition_test
        ;;
    A0_UNIT)
        test_unit a_operator_test
        ;;
    B0_UNIT)
        test_unit b_operator_test
        ;;
    C0_UNIT)
        test_unit c_operator_test
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
    ALL)
        # TODO run all tests
        ;;
    *)
        echo 'did not specify test...'
        ;;
esac
