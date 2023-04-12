#!/bin/bash

APIKEY=$YOUR_API_KEY

while read f;                    
do
name=$(echo $f | cut -d ',' -f 1)
ncbi-acc-download -F genbank $name --api-key $APIKEY 
biosample=$(grep "BioSample"  "$f".gbk | cut -d ':' -f 2 | tr -d ' ')
mv "$f".gbk gbk
echo $f,$biosample >> id_and_biosample.csv
echo $f
done < ../accessions.txt
