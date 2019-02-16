cd $MOR_DIR/cases/baf
$MOR_DIR/bin/linkm

echo 'baf'     > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME

$MOR_DIR/bin/gsnaps baf

if [ ${TEST: -4} == "UNIT" ]; then
    bash $MOR_DIR/tests/unit.sh $TEST
else
    case "$TEST" in
        BAF_INTEG)
            $MOR_DIR/tests/baf-test
            ;;
    esac
fi
