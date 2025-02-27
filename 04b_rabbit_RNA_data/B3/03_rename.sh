#!/bin/bash

for entry in *.fastq.gz
do
 mv $entry B3_$entry
done