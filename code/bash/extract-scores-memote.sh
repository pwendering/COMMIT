#!/bin/bash
# Extract memote scores from html report files.
# Prints the result to stdout.

get_score() {
  # Input 1 is filename
  # Input 2 is query
  echo $(cat $1 | grep -A2 "$2" | grep -o "[0-9][0-9\.\,]*")
}

# isolate ID
ID=$(head -1 $1)

# stoichiometric consistency
s_con=$(get_score $1 "Stoichiometric Consistency")

# mass balance
m_bal=$(get_score $1 "Mass Balance")

# charge balance
c_bal=$(get_score $1 "Charge Balance")

# metabolite connectivity
m_conn=$(get_score $1 "Metabolite Connectivity")

# unbounded flux
ubd_flx=$(get_score $1 "Unbounded Flux In Default Medium")

# stoichiometrically balanced cycles
sbc=$(get_score $1 "Stoichiometrically Balanced Cycles")

# universally blocked reactions
blocked=$(get_score $1 "Universally Blocked Reactions")

# consistency	subscore
cons_sub=$(get_score $1 "Sub Total" | grep -m 1 -o "[0-9][0-9\.]*" | head -n 1)

# total score
tot_score=$(get_score $1 "Total Score" | grep -m 1 -o "[0-9][0-9\.]*" | head -n 1)


echo -e "$ID\t$s_con\t$m_bal\t$c_bal\t$m_conn\t$ubd_flx\t$sbc\t$blocked\t$cons_sub\t$tot_score"
