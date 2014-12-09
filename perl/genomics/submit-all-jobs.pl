#!/bin/sh
ls *.q | perl -ne 'print "qsub $_"' | /bin/sh >& /dev/null


