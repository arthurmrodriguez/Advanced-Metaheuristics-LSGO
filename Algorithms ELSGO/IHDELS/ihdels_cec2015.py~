#!/usr/bin/python
#     pylint: disable=E1101
import argparse
import sys

from os import path

from ea import DE
from ea.PoolProb import PoolProb
from ea import DEcrossover
from ea.DE import EAresult

from scipy.optimize import fmin_l_bfgs_b, minimize

from numpy.random import seed, permutation, uniform, randint

from numpy import repeat, tile, ones, zeros, arange, array, array_str, clip, copy

from os.path import isfile

from functools import partial

from mts import mtsls

def get_improvement(alg_name, before, after):
    """
    Print the improvement with an algorithm
    """
    return "{0}: {1:.3e} -> {2:.3e} [{3:.2e}]\n".format(alg_name, before, after, before-after)

SR_global_MTS = []
SR_MTS = []

def apply_minimize(name, method, fitness_fun, bounds, current_best, current_best_fitness, maxevals, fid):
    global SR_MTS
    global SR_global_MTS

    lower = bounds[0][0]
    upper = bounds[0][1]

    if method == 'grad':
        sol, fit, info = fmin_l_bfgs_b(fitness_fun, x0=current_best, approx_grad=1, bounds=bounds, maxfun=maxevals, disp=False)
        funcalls = info['funcalls']
    elif method == 'powell':
        res = minimize(fitness_fun, x0=current_best, method='Powell', options={'disp': False, 'maxfev': maxevals})
        sol = clip(res.x, lower, upper)
        fit = fitness_fun(sol)
        funcalls = maxevals
    elif method == 'mts':
#        import ipdb
#        ipdb.set_trace()
        if name.lower() == "global":
            SR = SR_global_MTS
        else:
            SR = SR_MTS

        res, SR_MTS = mtsls(fitness_fun, current_best, current_best_fitness, lower, upper, maxevals, SR)
        sol = res.solution
        fit = res.fitness
        funcalls = maxevals
    else:
        raise NotImplementedError(method)

    if fit <= current_best_fitness:
        fid.write(get_improvement("{0} {1}".format(method.upper(), name), current_best_fitness, fit))
        return EAresult(solution=array(sol), fitness=fit, evaluations=funcalls)
    else:
        return EAresult(solution=current_best, fitness=current_best_fitness, evaluations=funcalls)

def random_population(lower, upper, dimension, size):
    return uniform(lower, upper, dimension*size).reshape((size, dimension))

def _group_fitness_eval(fitness, sol, group, values):
    """
    Eval the fitness of a solution changing the group variables for the values. The idea is to uses it with
    functools.partial to give each EA component a different fitness function

    :param fitness: function to apply.
    :param sol: reference solution.
    :param group: group of variables to evaluate.
    :param values: values assigned to the group of variables
    """
    sol[group] = values
    return fitness(sol)

def get_fitness_partial(fitness_fun, current_best_solution, variables):
    return partial(_group_fitness_eval, fitness_fun, current_best_solution, variables)

def get_group_variables(dim, partial_dim, blocksize):
    """
    Return a group of variables, random
    """
    size = partial_dim/blocksize
    num_blocks = dim/blocksize
    blocks = permutation(num_blocks).reshape((num_blocks/size, size))
    extend = repeat(blocks, blocksize, axis=1)*blocksize
    added = tile(arange(blocksize), size)
    return extend+added

def applySADE(crossover, fitness, funinfo, partial_dim, evals, population, current_best, fid):
    result = DE.DE(run_info=funinfo, replace=False, debug=False, crossoverFunction=crossover, dimension=partial_dim, F='a', CR='n', name_output=None,
            population=population, fun=fitness, max_evals=evals, run=1, initial_solution=current_best.solution)
    fid.write(get_improvement("DE partial", current_best.fitness, result.fitness))
    return result

optimo = True

def check_evals(totalevals, evals, bestFitness, globalBestFitness, fid):
    if not evals:
        return evals
    elif totalevals >= evals[0]:
        best = min(bestFitness, globalBestFitness)
        fid.write("[%.1e]: %e,%d\n" %(evals[0], best, totalevals))
        fid.flush()
        evals.pop(0)

    return evals

def debfgs(fitness_fun, funinfo, dim, evals, fid, 
        blocksize=10, partial_dim=20, evals_ls=500, evals_de=1000, evals_gs=5000, debug=False, adapted=False):
    """
    Implementation of the proposal for CEC2015
    """
    lower = funinfo['lower']
    upper = funinfo['upper']
    evals = evals[:]

    initial_sol = ones(dim)*((lower+upper)/2.0)
    current_best_fitness = fitness_fun(initial_sol)

    maxevals = evals[-1]
    totalevals = 1

    group_variables = get_group_variables(dim, partial_dim, blocksize)

    bounds = zip(ones(dim)*lower, ones(dim)*upper)
    bounds_partial = zip(ones(partial_dim)*lower, ones(partial_dim)*upper)

    current_best = EAresult(solution=initial_sol, fitness=current_best_fitness, evaluations=totalevals)
    popsize = min(partial_dim, 100)
    population = random_population(lower, upper, dim, popsize)
    crossover = DEcrossover.SADECrossover(2)
    best_global_solution = current_best.solution 
    best_global_fitness = current_best.fitness
    current_best_solution = best_global_solution

    max_de_stucked = 5
    applied_restart = False
    de_stucked = 0

    apply_de = apply_ls = True
    global SR_MTS

    if not adapted:
        method = 'grad'
    elif adapted == 1:
        pool = PoolProb(['de', 'grad'])
    elif adapted == 2:
        pool = PoolProb(['grad', 'powell'])
    
    global SR_global_MTS
    SR_global_MTS = ones(dim)*(upper-lower)*0.2
    pool = PoolProb(['mts', 'grad'])
    pool_global = PoolProb(['mts', 'grad'])

    while totalevals < maxevals:
        group_variables = get_group_variables(dim, partial_dim, blocksize)
        evals_gs = max(50*dim, 25000)
        evals_de = max(50*partial_dim, 25000)
        evals_ls = max(10*partial_dim, 5000)

        if not pool_global.is_empty():
            previous = current_best.fitness
            method_global = pool_global.get_new()
            current_best = apply_minimize("Global", method_global, fitness_fun, bounds, current_best_solution, current_best.fitness, evals_gs, fid)
            totalevals += current_best.evaluations
            pool_global.improvement(method_global, previous - current_best.fitness, 2)
            evals = check_evals(totalevals, evals, current_best.fitness, best_global_fitness, fid)
            current_best_solution = current_best.solution
            current_best_fitness = current_best.fitness

        for i, variables in enumerate(group_variables):
             best_global_fitness = fitness_fun(best_global_solution)
             fitness_partial = get_fitness_partial(fitness_fun, current_best_solution, variables)
             current_best = EAresult(solution=current_best_solution[variables], fitness=current_best_fitness, evaluations=0)
             previous_fitness = current_best_fitness
             SR_MTS = SR_global_MTS[variables]

             if adapted:
                method = pool.get_new()

                if adapted == 1:
                    apply_de = (method == 'de')
                    apply_ls = (method != 'de')
                
             if apply_de:
                result = applySADE(crossover, fitness_partial, funinfo, partial_dim, evals_de, population[:,variables], current_best, fid)
                improvement = current_best.fitness - result.fitness
                totalevals += result.evaluations
                evals = check_evals(totalevals, evals, result.fitness, best_global_fitness, fid)
                current_best = result

                if improvement:
                    de_stucked = 0
                else:
                    de_stucked += 1

                    if de_stucked >= max_de_stucked:
                        de_stucked = 0

                        if applied_restart:
                            apply_de = False
                        else:
                            population[:,variables] = random_population(lower, upper, partial_dim, popsize)
                            applied_restart = True

             if apply_ls:
                result = apply_minimize("Local", method, fitness_partial, bounds_partial, current_best.solution, current_best.fitness, evals_ls, fid)
                improvement = current_best.fitness - result.fitness
                totalevals += result.evaluations
                evals = check_evals(totalevals, evals, result.fitness, best_global_fitness, fid)
                current_best = result

             if adapted:
                pool.improvement(method, improvement, freq_update=10, minimum=.25)
#                fid.write("Prob[de,grad] : {0}\n".format(array_str(pool.get_prob(), precision=2)))


             # Restart if it is not improved
             if previous_fitness == result.fitness:
                random_ind = population[randint(popsize)]
                current_best = EAresult(solution=random_ind[variables], fitness=fitness_fun(random_ind), evaluations=0)

             current_best_solution[variables] = current_best.solution
             current_best_fitness = current_best.fitness
    
             if current_best_fitness < best_global_fitness:
                 best_global_fitness = current_best_fitness
                 best_global_solution = copy(current_best_solution)

             fid.write("{0:.2e}: with {1:d} evaluations\n".format(current_best.fitness, totalevals))
#             fid.write("improvement_group[{}] : {:.2e}\n".format(i, (initial_fitness - result.fitness)))
             fid.flush()

             if totalevals >= maxevals:
                break

    fid.write("%e,%s,%d\n" %(abs(best_global_fitness), ' '.join(map(str, best_global_solution)), totalevals))
    fid.flush()
    return result

from cec2013lsgo.cec2013 import Benchmark

def main(args):
    """
    Main program. It uses
    """
    parser = argparse.ArgumentParser(description=
"""
Run DE for experiments. F, CR must be float, or 'n' as a normal
""")
    parser.add_argument("-f", required=True, type=int, choices=range(1, 16), dest="function", help='function')
    parser.add_argument("-v", default=False, dest="verbose", action='store_true', help='verbose mode')
    parser.add_argument("-o", default=False, dest="override", action='store_true', help='override mode')
    parser.add_argument("-s", required=True, type=int, dest="seed", choices=range(1, 6), help='seed (1 - 5)')
    parser.add_argument("-r", default=5, type=int, dest="run", help='runs')
    parser.add_argument("-e", required=False, type=int, dest="maxevals", help='maxevals')
    parser.add_argument("-g", default=500, type=int, dest="groupsize", help='group_size')
    parser.add_argument("-a", default=3, type=int, dest="adapted", help='adapted', choices=range(4))

    #seeds
    seeds = [23, 45689, 97232447, 96793335, 12345679]

    args = parser.parse_args(args)
    fun = args.function
    dim = 1000

    if (args.maxevals):
        evals = map(int, [1.2e5, 6e5, 3e6])[:args.maxevals]
    else:
        evals = map(int, [1.2e5, 6e5, 3e6])

    bench = Benchmark()
    maxfuns = bench.get_num_functions()
    funinfo = bench.get_info(fun)

    if not (1 <= fun <= maxfuns and 1 <= args.seed <= 5):
        parser.print_help()
        sys.exit(1)

    if not (args.groupsize >= 10 and (dim % args.groupsize)==0):
        parser.print_help()
        sys.exit(1)

    if args.adapted:
        name = "L{0}DEBFGS".format(args.adapted)
    else:
        name = "DEBGFS"

    fname = name +"_reg{args.groupsize}_F{args.function}_{args.seed}r{args.run}.txt".format(args=args);

    output = path.join("results", fname)

    if not args.verbose and isfile(output) and not args.override:
        fin = open(output, 'rb')
        lines = fin.readlines()
        fin.close()

        if lines:
            return

    if not args.verbose:
        fid = open(output, 'w')
    else:
        fid = sys.stdout

    # Parameter commons
    fitness_fun = bench.get_function(fun)

    seed(seeds[args.seed-1])

    for _ in range(args.run):
        debfgs(fitness_fun, funinfo, dim, evals, fid, adapted=args.adapted, partial_dim=args.groupsize)

    fid.close()

if __name__ == '__main__':
    main(sys.argv[1:])
