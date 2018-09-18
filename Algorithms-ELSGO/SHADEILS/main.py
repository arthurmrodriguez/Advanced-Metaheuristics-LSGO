import shadeils as SHADEILS
import argparse
import sys
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
    description = "IHDELS"
    parser = argparse.ArgumentParser(description)
    parser.add_argument("-f", required=True, type=int, choices=range(1, 16), dest="function", help='function')
    parser.add_argument("-v", default=False, dest="verbose", action='store_true', help='verbose mode')
    parser.add_argument("-s", default=1, type=int, dest="seed", choices=range(1, 6), help='seed (1 - 5)')
    parser.add_argument("-r", default=5, type=int, dest="run", help='runs')
    parser.add_argument("-e", required=False, type=int, dest="maxevals", help='maxevals')
    parser.add_argument("-t", default=0.01, type=float, dest="threshold", help='threshold')
    parser.add_argument("-p", default=100, type=int, dest="popsize", help='population size')
    parser.add_argument("-H", default=None, type=int, dest="shade_h", help='SHADE history size')
    parser.add_argument("-d", default="results", type=str, dest="dir_output", help='directory output')
    #seeds
    seeds = [23, 45689, 97232447, 96793335, 12345679]
    args = parser.parse_args(args)
    fun = args.function
    dim = 1000

    print("Function: {0}".format(fun))
    print("Seed: {0}".format(args.seed))
    print("Treshold: {0}".format(args.threshold))
    print("Popsize: {0}".format(args.popsize))

    if args.shade_h is None:
        args.shade_h = min(args.popsize, 100)

    print("SHADE_H: {0}".format(args.shade_h))

    if (args.maxevals):
        evals = list(map(int, [1.2e5, 6e5, 3e6])[:args.maxevals])
    else:
        evals = list(map(int, [1.2e5, 6e5, 3e6]))

    print(evals)
    bench = Benchmark()
    maxfuns = bench.get_num_functions()
    funinfo = bench.get_info(fun)

    if not (1 <= fun <= maxfuns and 1 <= args.seed <= 5):
        parser.print_help()
        sys.exit(1)

    name = "SHADEILS"

    fname = name +"_pop{args.popsize}_H{args.shade_h}_t{args.threshold:.2f}_F{args.function}_{args.seed}r{args.run}.txt".format(args=args);

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
    fitness_fun = bench.get_function(fun)

    seed(seeds[args.seed-1])

    for _ in range(args.run):
        SR_MTS = []
        SR_global_MTS = []
        SHADEILS.ihshadels(fitness_fun, funinfo, dim, evals, fid, threshold=args.threshold, popsize=args.popsize, info_de=args.shade_h)
        bench.next_run()

    fid.close()

if __name__ == '__main__':
    main(sys.argv[1:])
