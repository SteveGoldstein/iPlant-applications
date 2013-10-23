User Guide
paml_run.pl

paml_run is a wrapper for Zieheng Yang's suite of Maximum Likelihood programs for phylogenetic analyses PAML. The wrapper is required to overcome incompatibilities between iPlant's execution and Discovery Environments and the way PAML handles options.
IN particular PAML control files codeml.ctl contains the paths for input and additional data files (e.g. replacement matrices). 
Because the programs are executed on the Condor execution node, these path are invalid and would result in the program failure. To overcome this, absolute paths are corrected so that they point at the correct files.
Furthermore, this wrapper script facilitate the building of user interfaces to modify the run parameters. These are specified in the wrapper command line with the syntax parameter=value.
Two parameters are always required: program and ctrlfile which specify which binary to run and which control file to use as a template. In most cases, the input seqfile and treefile also need to be specified.
To run codeml the proper syntax is
./paml_run.pl program=codeml ctrlfile=/User/path/codeml.ctl seqfile=/User/data/path/seqs.fas treefile=/User/data/path/tree.phy

This will generate a new control file based on /User/path/codeml.ctl but with the corrected file paths and then will run codeml using that file as input.