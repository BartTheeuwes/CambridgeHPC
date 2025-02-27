#!/bin/bash

for entry in *.fastq.gz
do
 mv $entry B4_$entry
done