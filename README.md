Inferring Genetic Architecture From Sibling Data 

Requirements: 2-column comma or whitespace separated text file 
Examples: 

Within Directory, Type: 

./sibArc.py data/trait1-test.csv --savePlot --out denovoExample
./sibArc.py data/trait2-test.csv --savePlot --out mendelianExample 


Outputs Produced: 


denovoExample.h2.out, mendelianExample.h2.out:            Sibling Inferred Strata Specific Heritability 
denovoExample.result.out, mendelianExample.result.out:    Inferred Genetic Architecture in Each Tail 
denovoExample.fig.png, mendelianExample.fig.png:          Genetic Architecture Plots 
