#!/bin/bash

set -v

make debug
bin/qvz -c 1 -f 0.5 -s codebook.txt test.in test.q
bin/qvz -x codebook.txt test.q test.dec
diff fref.txt test.dec
