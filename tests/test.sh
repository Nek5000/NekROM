# MOR Tests

cd $ROOT_DIR/cases/baf
$ROOT_DIR/bin/linkc

echo 'baf'     > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME

$ROOT_DIR/bin/gnaps baf

case "$TEST" in
    GRAMMIAN_UNIT)
        cp $ROOT_DIR/tests/grammian_test.f t.f
        ./nek5000
        ;;
    BAFF_COMP)
        $ROOT_DIR/tests/comp-baff-test
        ;;
    BAF_COMP)
        $ROOT_DIR/tests/comp-baf-test
        ;;
    *)
        echo 'did not specify test...'
        ;;
esac
