#!/bin/bash

set -e
cd fastqc_output_all_2/

echo "Unzipping..."
for filename in *.zip
    do
    unzip $filename
    done

echo "Saving summary..."
cat */summary.txt > fastqc_summaries.txt