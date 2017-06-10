#!/bin/bash

echo "      date     time $(free -m | grep total | sed -E 's/^    (.*)/\1/g')" > memory-usage.log
while true; do
    echo "$(date '+%Y-%m-%d %H:%M:%S') $(free -m | grep Mem: | sed 's/Mem://g')" >> memory-usage.log
    sleep 1
done
