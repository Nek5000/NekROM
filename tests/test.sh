# MOR Tests

test_unit() {
test=$(echo $1 | perl -ne 'print lc')
$ROOT_DIR/tests/test_template.sh $test
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

if [ ${TEST: -4} == "UNIT" ]; then
    test_unit $TEST
else
    case "$TEST" in
        BAF_INTEG)
            $ROOT_DIR/tests/baf-test
            ;;
    esac
fi
