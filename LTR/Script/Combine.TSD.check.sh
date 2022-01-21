#!/usr/bin/bash
for i in *.dir
do 
	cat $i/perfect.TSD.list >>perfect.TSD.list;
	cat $i/1mismatch.TSD.list >>1mismatch.TSD.list;
	cat $i/No_TSD_found.list >>No_TSD_found.list;
done
