#!/bin/bash
# Reproduce Gomila's first 2,000 finite windows (N = 690988 .. 692987) with his exact shard boundaries and arguments.
cd "$(dirname "$0")"
run() { ./build/triangle "$1" "$2" 16125 16125 100000 350708 10000000 6 220 14 0.005 t > "repro/p235711_$1_$2.log" 2>&1; }
for s in "690988 690988" "690989 690990" "690991 690995" "690996 691010" "691011 691050" "691051 691150" "691151 691500" "691501 693000"; do
  run $s &
done
wait
echo REPRO_DONE
