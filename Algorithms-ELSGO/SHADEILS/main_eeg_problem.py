import shadeils as SHADEILS
import argparse
import sys
import time
from os import path
from os.path import isfile
from cec2013lsgo.cec2013 import Benchmark
from numpy.random import seed


def main(args):
    global SR_MTS, SR_global_MTS

    """
    Main program. It uses
    Run DE for experiments. F, CR must be float, or 'n' as a normal
    """
    description = "SHADEILS"
    parser = argparse.ArgumentParser(description)
    parser.add_argument("-f", default=16, type=int, dest="function", help='function')
    # Argument for problem type: D4, D12 or D19, both with or without noise
    parser.add_argument("-pr",required=True, type=str, choices=list(("D4","D12","D19","D4N","D12N","D19N")),dest="problem",help='problem type')
    parser.add_argument("-v", default=False, dest="verbose", action='store_true', help='verbose mode')
    parser.add_argument("-s", default=1, type=int, dest="seed", choices=range(1, 6), help='seed (1 - 5)')
    parser.add_argument("-r", default=1, type=int, dest="run", help='runs')
    parser.add_argument("-e", required=False, type=int, dest="maxevals" ,help='maxevals')
    # Argument for index of evals, that is [1,2,3] -> [1.2e5, 6e5,3e6]
    parser.add_argument("-i", default=1, type=int, dest="index",choices=range(1,4),help="index of evals")
    parser.add_argument("-t", default=0.001, type=float, dest="threshold", help='threshold')
    parser.add_argument("-p", default=100, type=int, dest="popsize", help='population size')
    parser.add_argument("-H", default=3, type=int, dest="shade_h", help='SHADE history size')
    parser.add_argument("-d", default="results", type=str, dest="dir_output", help='directory output')

    # Time stamp for comparative study
    time0 = time.clock()

    #seeds
    seeds = [24, 45689, 97232447, 96793335, 12345679]
    args = parser.parse_args(args)
    fun = args.function
    problem = args.problem

    # Noise problem, False by default
    noise = False

    if(problem in ["D4","D4N"]):
        dim = 1024
    elif(problem in ["D12","D12N"]):
        dim = 3072
    elif(problem in ["D19","D19N"]):
        dim = 4864

    if("N" in problem):
        print("Problem with noise")
        noise = True

    print("EEG Problem: {0}".format(problem))
    print("Problem: {0}".format(problem))
    print("Dimension: {0}".format(dim))
    print("Seed: {0}".format(args.seed))
    print("Treshold: {0}".format(args.threshold))
    print("Popsize: {0}".format(args.popsize))

    if args.shade_h is None:
        args.shade_h = min(args.popsize, 100)

    print("SHADE_H: {0}".format(args.shade_h))

    if (args.maxevals):
        evals = list(map(int, [1.0e5, 6e5, 3e6])[:args.maxevals])
    else:
        evals = list(map(int, [1.0e5, 6e5, 3e6]))

    evals_index = args.index

    bench = Benchmark()
    maxfuns = bench.get_num_functions()
    funinfo = bench.get_info(fun,dim)

    if not (1 <= fun <= maxfuns and 1 <= args.seed <= 5):
        parser.print_help()
        sys.exit(1)

    name = "SHADEILS"

    fname = name + "_EEGProblem_" + problem + "_" + str(evals[evals_index-1]) +".txt"

    output = path.join(args.dir_output, fname)

    if not args.verbose and isfile(output):
        fin = open(output, 'rb')
        lines = fin.readlines()
        fin.close()

        if lines:
            print("Experiment already exists. Check results folder")
            return

    if not args.verbose:
        fid = open(output, 'w')
    else:
        fid = sys.stdout


    # Parameter commons
    #bench.set_algname("shadeils_restart0.1_pos")
    fitness_fun = bench.get_function(fun,dim,noise)

    seed(seeds[args.seed-1])

    for _ in range(args.run):
        SR_MTS = []
        SR_global_MTS = []
        SHADEILS.ihshadels(fitness_fun, funinfo, dim, evals,evals_index, fid, time0, threshold=args.threshold, popsize=args.popsize, info_de=args.shade_h)
        bench.next_run()

    fid.close()

if __name__ == '__main__':
    main(sys.argv[1:])
