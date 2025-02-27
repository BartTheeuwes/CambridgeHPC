#!/bin/bash

for entry in *.fastq.gz
do
 mv $entry B2_$entry
done