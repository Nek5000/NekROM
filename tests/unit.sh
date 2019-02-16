if [[ $1 =~ _L2_ ]]; then
    ifl2='.TRUE.'
    $MOR_DIR/bin/gops baf-l2
    wget https://uofi.box.com/shared/static/4u3ez373wzsmcv2exfw1dya2ex48152i.gz -O baf-l2-bas.tar.gz
    tar -xvzf baf-l2-bas.tar.gz
    mv baf-l2-bas/* .
else
    ifl2='.FALSE.'
    $MOR_DIR/bin/gops baf-h10
    wget https://uofi.box.com/shared/static/1br5da8zf9ldqpwut477589621vzrxw6.gz -O baf-h10-bas.tar.gz
    tar -xvzf baf-h10-bas.tar.gz
    mv baf-h10-bas/* .
fi

mkdir ops

name="$(echo $1 | perl -pe 's/_(L2|H10)_/_/g')(${ifl2})"
$MOR_DIR/tests/test_template.sh $(echo $name | perl -ne 'print lc')

ls ../../data/baf/baf0.f* > file.list
ls ../../data/baf/avgbaf0.f* > avg.list

ln -s $MOR_DIR/data/baf/baf0.f00001 r0.f00001

$SOURCE_ROOT/bin/makenek test

$SOURCE_ROOT/bin/genmap << Z
baf
.01
Z

./nek5000
