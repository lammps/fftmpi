#!/bin/csh
# generate a PDF version of Manual

txt2html -b *.txt

htmldoc --title --toctitle "Table of Contents" --tocfooter ..i --toclevels 4 --header ... --footer ..1 --size letter --linkstyle plain --linkcolor blue -f Manual.pdf Manual.html quick.html intro.html compile.html buildtest.html runtest.html layout.html usage.html api.html api_create.html api_setup.html api_tune.html api_compute.html api_stats.html api_remap.html

txt2html *.txt
