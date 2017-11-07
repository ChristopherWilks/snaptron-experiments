#!/usr/bin/env sh

sh ./examples.sh > examples.sh.test 2>&1
diff examples.sh.test tests/examples.sh.output
