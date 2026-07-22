#!/bin/bash

> resultats_tous.txt

for cas in 1 2 3
do
  for h in 2 4 8 16 32 64
  do
    echo "===== cas=$cas, mesh=d1q1_$h =====" >> resultats_tous.txt
    echo -e "$cas\nd1q1_$h\nq\nq\nq" | ./prog >> resultats_tous.txt 2>&1
    echo "" >> resultats_tous.txt
  done
done

for cas in 4 5 6
do
  for h in 2 4 8 16 32 64
  do
    echo "===== cas=$cas, mesh=d2q1_$h =====" >> resultats_tous.txt
    echo -e "$cas\nd2q1_$h\nq\nq\nq" | ./prog >> resultats_tous.txt 2>&1
    echo "" >> resultats_tous.txt
  done
done
