#!/bin/csh
# generate a PDF version of Manual

txt2html -b *.txt

htmldoc --title --toctitle "Table of Contents" --tocfooter ..i --toclevels 4 --header ... --footer ..1 --size letter --linkstyle plain --linkcolor blue -f Manual.pdf Manual.html quick.html intro.html features.html model.html usage.html launch.html buildlib.html testapps.html buildtest.html runtest.html api.html language.html create.html exchange.html query.html python.html errors.html

txt2html *.txt
