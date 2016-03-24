## Installation

1. Follow the instruction [1][1] and [2][2] to install EDG 4.x-based ROSE.
2. Edit `set.rose` to set the correct PATH to your install and run `. set.rose`. 
3. `$ cd tools/onto-build-rose/` and run `make`

[1]: https://en.wikibooks.org/wiki/ROSE_Compiler_Framework/Installation
[2]: https://en.wikibooks.org/wiki/ROSE_Compiler_Framework/Virtual_Machine_Image#V2

## Directory

+ ontology/ - owl files.
+ projects/ - contains separate analysis projects.
+ test/ - test bench and run test script.
+ tools/ - the knowledge generator: use ROSE frontend to parse the C code and build knowledge base.

## Usage and detail

### Generate ontology represented database from program

	$ rosePrgKnowledgeBuilder.exe -c -w -emit-owl out.owl input.c [-I/extra-include-dir]


### To run a prolog inference

Canonical Loop

	$ swipl --nosignal --quiet projects/canonicalloop/run.pl out.owl report.txt

or CFG

	$ swipl --nosignal --quiet projects/cfg_test/run.pl out.owl report.txt
