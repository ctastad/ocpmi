#!/bin/bash

source full_send_vars.txt


parallel 'Rscript --vanilla ../dnacopy_seg_plot.R {}' ::: *.copynumber
