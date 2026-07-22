#!/bin/bash

> erreurs.txt

for cas in 1 2 3
do
  for h in 2 4 8 16 32 64
  do
    echo "===== cas=$cas d1t1_$h =====" >> erreurs.txt
    printf "%s\n%s\nq\n" "$cas" "d1t1_$h" | ./prog 2>&1 | grep -i "erreur" >> erreurs.txt
    echo "" >> erreurs.txt
  done
done

for cas in 4 5 6
do
  for h in 2 4 8 16 32 64
  do
    echo "===== cas=$cas d2t1_$h =====" >> erreurs.txt
    printf "%s\n%s\nq\n" "$cas" "d2t1_$h" | ./prog 2>&1 | grep -i "erreur" >> erreurs.txt
    echo "" >> erreurs.txt
  done
done
