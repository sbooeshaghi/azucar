#!/usr/bin/env bash
gunzip out/*.gz
../../assign/mx_index.py assign/markers.txt assign/
../../assign/mx_select.py assign/markers.txt out/genes.txt assign/
../../assign/mx_extract.py out/matrix.mtx assign/select.txt assign/matrix_select.mtx
../../assign/mx_assign.py assign/matrix_select.mtx assign/markers.ec assign/assignments.txt.gz out/barcodes.txt assign/groups.txt
gzip out/dbco.txt out/matrix.mtx out/genes.txt out/barcodes.txt
