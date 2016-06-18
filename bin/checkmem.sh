ps -orss= -u `whoami` | awk '{s+=$1}END{print s}'
