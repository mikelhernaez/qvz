#!/bin/bash

set -v

make debug
bin/qvz -c 1 -f 0.5 -s test.in test.q
bin/qvz -x test.q test.dec
diff fref.txt test.dec
