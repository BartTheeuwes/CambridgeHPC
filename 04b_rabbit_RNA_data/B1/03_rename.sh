#!/bin/bash

for entry in *.fastq.gz
do
 mv $entry B1_$entry
done