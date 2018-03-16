===================
About this document
===================

This project contains the Source code used for the work:

"Molina, Daniel; Herrera, Francisco, "Iterative hybridization of DE with local search
for the CEC'2015 special session on large scale global optimization," in Evolutionary
Computation (CEC), 2015 IEEE Congress on , vol., no., pp.1974-1978, 25-28 May 2015 doi:
10.1109/CEC.2015.7257127"

LICENSE
-------
The software is free software, in particular it is under the GNU v3.0 LICENSE.

Requirements
------------

This software has been developed and tested with:

- Python 2.7.9. 
- numpy 1.8.2.
- Cython 0.21.7.

However, other versions of these packages could also be used. 

Install
-------

1- Compile and Install the package ea:

$ cd ea
$ python setup.py install --user
$ cd ..

2- Use the program ihdels_cec2015.py with:

$ python ihdels_cec2015.py ...
or
$ ./ihdels_cec2015.py ... (if 

Running
-------

The program run with ihdels_cec2015.py. All parameters are the same than the presented in the IEEE CEC'2015
congress. Only two parameters are required:

- ./ihdels_cec2015.py -f <fun> -s <seed>

Required parameters:

- <fun> is the function to evaluate (between [1-15]).

- <seed> is a number (in [1-5]) that represents the seed value used
  for the random functions.

Aditional parameters:

- r <run> is the number of total run.

Note: The program stores all its data in the directory results/, thus you should check that this directory exists
before running the program. 

Understanding results
----------------------
In results/ the program stores a file with the information about the run.
There is a lot of information during the run. If you are only interested in the final results, you can consider
only the lines with the format:

[evals]: best_value, evals_best_value

where evals implies the number of evaluations (1.2e5, 6.0e+5, 3.0e+6),
best_value is the best obtained fitness value, and evals_best_value is the evaluation in which
that best value was obtained.

Tip:
$ grep "^\[" filename
