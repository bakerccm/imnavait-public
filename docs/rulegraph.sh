#!/bin/bash

# snakemake -c 1 --use-conda --forceall --rulegraph | dot -Tpdf >docs/rulegraph.pdf
# snakemake -c 1 --use-conda --forceall --rulegraph | dot -Tpng >docs/rulegraph.png

snakemake -c 1 --use-conda --forceall --rulegraph >docs/temp
dot -Tpdf docs/temp >docs/rulegraph.pdf
dot -Tpng docs/temp >docs/rulegraph.png
rm docs/temp
