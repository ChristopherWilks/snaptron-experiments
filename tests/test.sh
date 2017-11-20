#!/usr/bin/env sh

python client/test_query_snaptron.py
sh ./examples.sh > examples.sh.test 2>&1
diff examples.sh.test tests/examples.sh.output
