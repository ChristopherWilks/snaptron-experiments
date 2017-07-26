#!/bin/bash

./scripts/run_ts.sh > run_ts.out
./scripts/run_jir.sh
./scripts/run_ssc.sh
./scripts/run_psi.sh
./scripts/run_intersection.sh
