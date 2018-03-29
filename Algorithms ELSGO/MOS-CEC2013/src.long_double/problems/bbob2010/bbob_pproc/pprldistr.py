#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""Creates run length distribution figures."""

from __future__ import absolute_import

import os
import numpy
import pickle
import matplotlib.pyplot as plt
from pdb import set_trace
from bbob_pproc import bootstrap
from bbob_pproc.ppfig import consecutiveNumbers

rldColors = ('k', 'c', 'm', 'r', 'k', 'c', 'm', 'r', 'k', 'c', 'm', 'r')
rldUnsuccColors = ('k', 'c', 'm', 'k', 'c', 'm', 'k', 'c', 'm', 'k', 'c', 'm')  # should not be too short

# Used as a global to store the largest xmax and align the FV ECD figures.
fmax = None
evalfmax = None
figformat = ('eps', 'pdf') # Controls the output when using the main method

filename = 'pprldistr2009_1e-8.pickle'
filename = os.path.join(os.path.split(__file__)[0], filename)
isAlgorithm2009Found = True
try:
    f = open(filename,'r')
    dict2009 = pickle.load(f)
except IOError, (errno, strerror):
    print "I/O error(%s): %s" % (errno, strerror)
    isAlgorithm2009Found = False
    print 'Could not find file: ', filename
else:
    f.close()

def plotECDF(x, n=None, plotArgs={}):
    if n is None:
        n = len(x)
    nx = len(x)
    if n == 0 or nx == 0:
        res = plt.plot([], [], **plotArgs)
    else:
        x2 = numpy.hstack(numpy.repeat(sorted(x), 2))
        y2 = numpy.hstack([0.0,
                           numpy.repeat(numpy.arange(1, nx) / float(n), 2),
                           float(nx)/n])
        res = plt.plot(x2, y2, **plotArgs)
    return res

def beautifyECDF(axish=None):
    if axish is None:
        axish = plt.gca()
    plt.ylim(0.0, 1.0)
    plt.yticks(numpy.array((0., 0.25, 0.5, 0.75, 1.0)),
               ('0.0', '', '0.5', '', '1.0'))
    axish.grid('True')

def beautifyRLD(figHandle, figureName, maxEvalsF, fileFormat=('pdf', 'eps'),
                text=None, verbose=True):
    """Format the figure of the run length distribution and save into files."""

    axisHandle = figHandle.gca()
    axisHandle.set_xscale('log')
    #plt.axvline(x=maxEvalsF, color='k')
    plt.xlim(1.0, maxEvalsF ** 1.05)
    axisHandle.set_xlabel('log10 of FEvals / DIM')
    axisHandle.set_ylabel('proportion of trials')
    # Grid options
    xtic = axisHandle.get_xticks()
    newxtic = []
    for j in xtic:
        newxtic.append('%d' % round(numpy.log10(j)))
    axisHandle.set_xticklabels(newxtic)

    beautifyECDF()

    plt.text(0.5, 0.93, text, horizontalalignment="center",
             transform=axisHandle.transAxes)
             #bbox=dict(ec='k', fill=False), 

    #set_trace()
    plt.legend(loc='best')
    #if legend:
        #axisHandle.legend(legend, locLegend)

    # Save figure
    for entry in fileFormat:
        plt.savefig(figureName + '.' + entry, dpi = 300,
                    format = entry)
        if verbose:
            print 'Wrote figure in %s.' %(figureName + '.' + entry)

def plotRLDistr(dsList, fvalueToReach, maxEvalsF, plotArgs={},
                 verbose=True):
    """Creates run length distributions from a sequence dataSetList.

    Keyword arguments:
    dataSetList
    fvalueToReach
    verbose

    Outputs:
    res -- resulting plot.
    fsolved -- number of different functions solved.
    funcs -- number of different function considered.
    """

    x = []
    nn = 0
    fsolved = set()
    funcs = set()
    for i in dsList:
        funcs.add(i.funcId)
        for j in i.evals:
            if j[0] <= fvalueToReach[i.funcId]:
                #set_trace()
                tmp = j[1:]
                x.extend(tmp[numpy.isfinite(tmp)]/float(i.dim))
                fsolved.add(i.funcId)
                #TODO: what if j[numpy.isfinite(j)] is empty
                break
        nn += i.nbRuns()
    #set_trace()
    kwargs = plotArgs.copy()
    try:
        label = ''
        if len(set(fvalueToReach.values())):
            label += '%+d:' % (numpy.log10(fvalueToReach[i.funcId]))
        label += '%d/%d' % (len(fsolved), len(funcs))
        kwargs['label'] = kwargs.setdefault('label', label)
    except TypeError: # fvalueToReach == 0. for instance...
        # no label
        pass

    #TODO: res = plotECDF(x, nn, kwargs) # Why not?
    n = len(x)
    if n == 0:
        res = plt.plot([], [], **kwargs)
    else:
        x.sort()
        x2 = numpy.hstack([numpy.repeat(x, 2), maxEvalsF ** 1.05])
        # maxEvalsF: used for the limit of the plot
        y2 = numpy.hstack([0.0,
                           numpy.repeat(numpy.arange(1, n+1)/float(nn), 2)])
        res = plt.plot(x2, y2, **kwargs)

    return res#, fsolved, funcs

def plotERTDistr(dsList, fvalueToReach, plotArgs=None, verbose=True):
    """Creates estimated run time distributions from a sequence dataSetList.

    Keyword arguments:
    dsList
    fvalueToReach
    verbose

    Outputs:
    res -- resulting plot.
    fsolved -- number of different functions solved.
    funcs -- number of different function considered.
    """

    x = []
    nn = 0
    samplesize = 1000 # samplesize is at least 1000
    percentiles = 0.5 # could be anything...

    for i in dsList:
        #funcs.add(i.funcId)
        for j in i.evals:
            if j[0] <= fvalueToReach[i.funcId]:
                runlengthsucc = j[1:][numpy.isfinite(j[1:])]
                runlengthunsucc = i.maxevals[numpy.isnan(j[1:])]
                tmp = bootstrap.drawSP(runlengthsucc, runlengthunsucc,
                                       percentiles=percentiles,
                                       samplesize=samplesize)
                x.extend(tmp[1])
                break
        nn += samplesize
    #set_trace()
    res = plotECDF(x, nn, plotArgs)

    return res

def generateRLData(evals, targets):
    """Determine the running lengths for attaining the targets.

    Keyword arguments:
    evals -- numpy array with the first column corresponding to the function
      values and the following columns being the number of function evaluations
      for reaching this function value
    targets -- target function values of interest

    Output:
    list of arrays containing the number of function evaluations for reaching
    the target function values in target.
    """

    res = {}
    it = reversed(evals) # expect evals to be sorted by decreasing function values
    prevline = numpy.array([-numpy.inf] + [numpy.nan] * (numpy.shape(evals)[1]-1))
    try:
        line = it.next()
    except StopIteration:
        # evals is an empty array
        return res

    for t in sorted(targets):
        while line[0] <= t:
            prevline = line
            try:
                line = it.next()
            except StopIteration:
                break
        #if prevline[0] > t:
            #prevline.copy()
        #set_trace()
        res[t] = prevline.copy() # is copy necessary?
    return res

def beautifyFVD(figHandle, figureName, fileFormat=('pdf', 'eps'),
                isStoringXMax=False, text=None, verbose=True):
    """Formats the figure of the run length distribution.

    Keyword arguments:
    isStoringMaxF -- if set to True, the first call BeautifyVD sets the global
                     fmax and all subsequent call will have the same maximum
                     xlim.
    """

    axisHandle = figHandle.gca()
    axisHandle.set_xscale('log')

    if isStoringXMax:
        global fmax
    else:
        fmax = None

    if not fmax:
        xmin, fmax = plt.xlim()
    plt.xlim(1., fmax)

    #axisHandle.invert_xaxis()
    axisHandle.set_xlabel('log10 of Df / Dftarget')
    # axisHandle.set_ylabel('proportion of successful trials')
    # Grid options
    beautifyECDF()

    xtic = axisHandle.get_xticks()
    newxtic = []
    for j in xtic:
        newxtic.append('%d' % round(numpy.log10(j)))
    axisHandle.set_xticklabels(newxtic)
    axisHandle.set_yticklabels(())

    plt.text(0.98, 0.02, text, horizontalalignment="right",
             transform=axisHandle.transAxes)
             #bbox=dict(ec='k', fill=False), 

    # Save figure
    for entry in fileFormat:
        plt.savefig(figureName + '.' + entry, dpi = 300,
                    format = entry)
        if verbose:
            print 'Wrote figure in %s.' %(figureName + '.' + entry)

def plotFVDistr(dataSetList, fvalueToReach, maxEvalsF, plotArgs={},
                 verbose=True):
    """Creates empirical cumulative distribution functions of final function
    values plot from a sequence of indexEntries.

    Keyword arguments:
    indexEntries -- sequence of IndexEntry to process.
    fvalueToReach -- float used for the lower limit of the plot
    maxEvalsF -- indicates which vertical data to display.
    verbose -- controls verbosity.

    Outputs: a plot of a run length distribution.
    """

    x = []
    nn = 0
    for i in dataSetList:
        for j in i.funvals:
            if j[0] >= maxEvalsF * i.dim:
                break

        tmp = j[1:].copy() / fvalueToReach[i.funcId]
        tmp[tmp<=0.] = 1.
        # TODO: HACK, is almost ok since the xmin in the figure is 1
        x.extend(tmp)
        nn += i.nbRuns()

    res = plotECDF(x, nn, plotArgs)

    return res

def comp(dsList0, dsList1, valuesOfInterest, isStoringXMax=False,
         outputdir='', info='default', verbose=True):
    """Generate figures of empirical cumulative distribution functions.
    Dashed lines will correspond to ALG0 and solid lines to ALG1.

    Keyword arguments:
    dsList0 -- list of DataSet instances for ALG0.
    dsList1 -- list of DataSet instances for ALG1
    valuesOfInterest -- target function values to be displayed.
    isStoringXMax -- if set to True, the first call BeautifyVD sets the globals
                     fmax and maxEvals and all subsequent calls will use these
                     values as rightmost xlim in the generated figures.
     -- if set to True, the first call BeautifyVD sets the global
                     fmax and all subsequent call will have the same maximum
                     xlim.
    outputdir -- output directory (must exist)
    info --- string suffix for output file names.

    Outputs:
    Image files of the empirical cumulative distribution functions.
    """

    plt.rc("axes", labelsize=20, titlesize=24)
    plt.rc("xtick", labelsize=20)
    plt.rc("ytick", labelsize=20)
    plt.rc("font", size=20)
    plt.rc("legend", fontsize=20)

    maxEvalsFactor = max(max(i.mMaxEvals()/i.dim for i in dsList0),
                         max(i.mMaxEvals()/i.dim for i in dsList1))

    if isStoringXMax:
        global evalfmax
    else:
        evalfmax = None

    if not evalfmax:
        evalfmax = maxEvalsFactor

    figureName = os.path.join(outputdir,'pprldistr2_%s' %(info))
    fig = plt.figure()
    legend = []
    for j in range(len(valuesOfInterest)):
        tmp = plotRLDistr(dsList0, valuesOfInterest[j], evalfmax,
                          verbose=verbose)

        if not tmp is None:
            plt.setp(tmp, 'color', rldColors[j])
            plt.setp(tmp, 'ls', '--')
            plt.setp(tmp, 'label', None) # Hack for the legend

            if rldColors[j] == 'r':  # 1e-8 in bold
                plt.setp(tmp, 'linewidth', 3)

        tmp = plotRLDistr(dsList1, valuesOfInterest[j], evalfmax,
                          verbose=verbose)

        if not tmp is None:
            plt.setp(tmp, 'color', rldColors[j])
            # Hack for the legend.
            plt.setp(tmp, 'label', ('%+d' % (numpy.log10(valuesOfInterest[j][1]))))

            if rldColors[j] == 'r':  # 1e-8 in bold
                plt.setp(tmp, 'linewidth', 3)

    funcs = set(i.funcId for i in dsList0) | set(i.funcId for i in dsList1)
    text = 'f%s' % (consecutiveNumbers(sorted(funcs)))

    if isAlgorithm2009Found:
        d = set.union(set(i.dim for i in dsList0),
                      set(i.dim for i in dsList1)).pop() # Get only one element...
        for alg in dict2009:
            x = []
            nn = 0
            try:
                tmp = dict2009[alg]
                for f in funcs:
                    tmp[f][d] # simply test that they exists
            except KeyError:
                continue

            for f in funcs:
                tmp2 = tmp[f][d][0][1:]
                # [0], because the maximum #evals is also recorded
                # [1:] because the target function value is recorded
                x.append(tmp2[numpy.isnan(tmp2) == False])
                nn += len(tmp2)

            if x:
                x.append([(evalfmax*d) ** 1.05])
                x = numpy.hstack(x)

                plotECDF(x[numpy.isfinite(x)]/d, nn,
                         {'color': 'wheat', 'ls': '-', 'zorder': -1})

    plt.axvline(max(i.mMaxEvals()/i.dim for i in dsList0), ls='--', color='k')
    plt.axvline(max(i.mMaxEvals()/i.dim for i in dsList1), color='k')
    beautifyRLD(fig, figureName, evalfmax, fileFormat=figformat, text=text,
                verbose=verbose)
    plt.close(fig)

    # The figures generated by the lines below are not displayed in the
    # BBOB2010 template: should they be done?
    #figureName = os.path.join(outputdir,'ppfvdistr2_%s' %(info))
    #fig = plt.figure()
    #for j in range(len(valuesOfInterest)):
        ##set_trace()
        #tmp = plotFVDistr(dsList0, valuesOfInterest[j],
                          #evalfmax, verbose=verbose)
        ##if not tmp is None:
        #plt.setp(tmp, 'color', rldColors[j])
        #plt.setp(tmp, 'ls', '--')
        #if rldColors [j] == 'r':  # 1e-8 in bold
            #plt.setp(tmp, 'linewidth', 3)

        #tmp = plotFVDistr(dsList1, valuesOfInterest[j],
                          #evalfmax, verbose=verbose)
        ##if not tmp is None:
        #plt.setp(tmp, 'color', rldColors[j])
        #plt.setp(tmp, 'ls', '--')
        #if rldColors [j] == 'r':  # 1e-8 in bold
            #plt.setp(tmp, 'linewidth', 3)

    #tmp = numpy.floor(numpy.log10(evalfmax))
    ## coloring left to right:
    ##maxEvalsF = numpy.power(10, numpy.arange(tmp, 0, -1) - 1)
    ## coloring right to left:
    #maxEvalsF = numpy.power(10, numpy.arange(0, tmp))

    ##The last index of valuesOfInterest is still used in this loop.
    ##set_trace()
    ## Plot lines for different number of function evaluations
    #for k in range(len(maxEvalsF)):
        #tmp = plotFVDistr(dsList0, valuesOfInterest[j],
                           #maxEvalsF=maxEvalsF[k], verbose=verbose)
        #plt.setp(tmp, 'color', rldUnsuccColors[k])
        #plt.setp(tmp, 'ls', '--')

        #tmp = plotFVDistr(dsList1, valuesOfInterest[j],
                           #maxEvalsF=maxEvalsF[k], verbose=verbose)
        #plt.setp(tmp, 'color', rldUnsuccColors[k])

    #beautifyFVD(fig, figureName, fileFormat=figformat, text=text,
                #isStoringXMax=isStoringXMax, verbose=verbose)

    #plt.close(fig)

    plt.rcdefaults()

def main(dsList, valuesOfInterest, isStoringXMax=False, outputdir='',
         info='default', verbose=True):
    """Generate figures of empirical cumulative distribution functions.

    Keyword arguments:
    dsList -- list of DataSet instances to process.
    valuesOfInterest -- target function values to be displayed.
    isStoringXMax -- if set to True, the first call BeautifyVD sets the globals
                     fmax and maxEvals and all subsequent calls will use these
                     values as rightmost xlim in the generated figures.
    outputdir -- output directory (must exist)
    info --- string suffix for output file names.

    Outputs:
    Image files of the empirical cumulative distribution functions.
    """

    plt.rc("axes", labelsize=20, titlesize=24)
    plt.rc("xtick", labelsize=20)
    plt.rc("ytick", labelsize=20)
    plt.rc("font", size=20)
    plt.rc("legend", fontsize=20)

    maxEvalsFactor = max(i.mMaxEvals()/i.dim for i in dsList)
    #maxEvalsFactorCeil = numpy.power(10,
                                     #numpy.ceil(numpy.log10(maxEvalsFactor)))

    if isStoringXMax:
        global evalfmax
    else:
        evalfmax = None

    if not evalfmax:
        evalfmax = maxEvalsFactor

    figureName = os.path.join(outputdir,'pprldistr_%s' %(info))
    fig = plt.figure()
    legend = []
    for j in range(len(valuesOfInterest)):
        tmp = plotRLDistr(dsList, valuesOfInterest[j], evalfmax,
                          verbose=verbose)

        if not tmp is None:
            plt.setp(tmp, 'color', rldColors[j])
            if rldColors[j] == 'r':  # 1e-8 in bold
                plt.setp(tmp, 'linewidth', 3)

    funcs = list(i.funcId for i in dsList)
    text = 'f%s' % (consecutiveNumbers(sorted(funcs)))

    if isAlgorithm2009Found:
        d = set(i.dim for i in dsList).pop() # Get only one element...
        for alg in dict2009:
            x = []
            nn = 0
            try:
                tmp = dict2009[alg]
                for f in funcs:
                    tmp[f][d] # simply test that they exists
            except KeyError:
                continue

            for f in funcs:
                tmp2 = tmp[f][d][0][1:]
                # [0], because the maximum #evals is also recorded
                # [1:] because the target function value is recorded
                x.append(tmp2[numpy.isnan(tmp2) == False])
                nn += len(tmp2)

            if x:
                x.append([(evalfmax*d) ** 1.05])
                x = numpy.hstack(x)

                plotECDF(x[numpy.isfinite(x)]/float(d), nn,
                         {'color': 'wheat', 'ls': '-', 'zorder': -1})

    plt.axvline(x=maxEvalsFactor, color='k')
    beautifyRLD(fig, figureName, evalfmax, fileFormat=figformat, text=text,
                verbose=verbose)

    plt.close(fig)

    figureName = os.path.join(outputdir,'ppfvdistr_%s' %(info))
    fig = plt.figure()
    for j in range(len(valuesOfInterest)):

        tmp = plotFVDistr(dsList, valuesOfInterest[j],
                          evalfmax, verbose=verbose)

        plt.setp(tmp, 'color', rldColors[j])
        if rldColors [j] == 'r':  # 1e-8 in bold
            plt.setp(tmp, 'linewidth', 3)

    tmp = numpy.floor(numpy.log10(evalfmax))
    # coloring left to right:
    #maxEvalsF = numpy.power(10, numpy.arange(tmp, 0, -1) - 1)
    # coloring right to left:
    maxEvalsF = numpy.power(10, numpy.arange(0, tmp))

    for k in range(len(maxEvalsF)):
        tmp = plotFVDistr(dsList, valuesOfInterest[-1],
                          maxEvalsF=maxEvalsF[k], verbose=verbose)
        plt.setp(tmp, 'color', rldUnsuccColors[k])

    beautifyFVD(fig, figureName, fileFormat=figformat, text=text,
                isStoringXMax=isStoringXMax, verbose=verbose)

    plt.close(fig)

    plt.rcdefaults()
