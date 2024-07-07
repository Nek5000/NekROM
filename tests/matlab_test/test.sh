cd $MOR_DIR/tests
source travis_fn
mkdir t
cd t
gdown https://drive.google.com/uc?id=1-5ryuIACyWtOzSl3y0y0SYKBqIGofrO3;;
unzip Re200.zip
cd Re200/T20/
cp $MOR_DIR/bin/matlab_driver/*.m .
cp $MOR_DIR/tests/matlab_test/test.m .
matlab -batch "run('test.m')"