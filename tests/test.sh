$MOR_DIR/bin/linkm

echo 'test'     > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME

$MOR_DIR/bin/gsnaps cyl_rect

if [ ${TEST: -4} == "UNIT" ]; then
    bash $MOR_DIR/tests/unit.sh $TEST
else
    case "$TEST" in
        BAF_INTEG)
            $MOR_DIR/tests/baf-test
            ;;
    esac
fi
