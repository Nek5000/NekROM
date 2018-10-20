# MOR Tests

cd $ROOT_DIR/cases/baf
$ROOT_DIR/bin/linkc

echo 'baf'     > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME

$ROOT_DIR/bin/gsnaps baf
$ROOT_DIR/bin/gops   baf

case "$TEST" in
#   GRAMMIAN_UNIT)
#       $ROOT_DIR/tests/grammian_test.sh
#       ;;
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
