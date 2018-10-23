if [[ $1 =~ _L2_ ]]; then
    ifl2='.TRUE.'
    $ROOT_DIR/bin/gops baf-l2
else
    ifl2='.FALSE.'
    $ROOT_DIR/bin/gops baf-h10
fi

name="$(echo $1 | perl -pe 's/_(L2|H10)_/_/g')(${ifl2})"
$ROOT_DIR/tests/test_template.sh $(echo $name | perl -ne 'print lc')

ls $ROOT_DIR/data/baf | sed 's/^/..\/..\/data\/baf\//g' > file.list

$SOURCE_ROOT/bin/makenek test

$SOURCE_ROOT/bin/genmap << Z
baf
.01
Z

./nek5000
