#!/bin/bash

for var in "$@"
do
	./calcTopas.sh "$var"
done

telegram-send --format markdown "*** ${HOSTNAME:11:7}: All queued simulations finished ***"
echo "All done!"
