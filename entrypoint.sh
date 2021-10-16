#!/bin/bash

export PYTHONPATH=/src

autoreconf -vif -I m4
./configure 
make install
make -j1 distcheck || exit 1

exec "$@"