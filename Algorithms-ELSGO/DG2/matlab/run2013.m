% Author: Mohammad Nabi Omidvar
% email address: mn.omidvar AT gmail.com
%
% ------------
% Description:
% ------------
% This files is the entry point for running the differential gropuing algorithm.

clear all;

% Specify the functions that you want the differential grouping algorithm to 
% identify its underlying grouping structure.

addpath('cec2013');
addpath('cec2013/datafiles');

% The function number for which a grouping is needed.
func_num = 9;

global initial_flag;
initial_flag = 0;

t1 = [1 4 7 8 11 12 13 14 15];
t2 = [2 5 9];
t3 = [3 6 10];

% setting the upper and lower bounds
if (ismember(func_num, t1))
    lb = -100;
    ub = 100;
elseif (ismember(func_num, t2))
    lb = -5;
    ub = 5;
elseif (ismember(func_num, t3))
    lb = -32;
    ub = 32;
end

opts.lbound  = lb;
opts.ubound  = ub;
opts.dim     = 1000;

if (ismember(func_num, [13 14]))
    opts.dim = 905;
end

% a wrapper function for the CEC LSGO suites.
objfun = @(x)(benchmark_func(x, func_num));

[delta, lambda, evaluations] = ism(objfun, opts);
[nonseps, seps, theta, epsilon] = dsm(evaluations, lambda, opts.dim);

filename = sprintf('./results/2013/F%02d.mat', func_num);
save (filename, 'delta', 'lambda', 'evaluations', 'nonseps', 'seps', 'theta', 'epsilon', '-v7');

