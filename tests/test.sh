cd $ROOT_DIR/cases/baf
$ROOT_DIR/bin/linkc

echo 'baf'     > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME

$ROOT_DIR/bin/gsnaps baf

if [ ${TEST: -4} == "UNIT" ]; then
    bash $ROOT_DIR/tests/unit.sh $TEST
else
    case "$TEST" in
        BAF_INTEG)
            $ROOT_DIR/tests/baf-test
            ;;
    esac
fi
