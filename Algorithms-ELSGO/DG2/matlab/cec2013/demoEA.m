clear;

tic;

NP = 100; % population size
runs = 2; % number of independent runs for each function, should be set to 25
Max_FEs = 1000; % maximal number of FEs, should be set to 3e+06
global initial_flag; % the global flag used in test suite
format short e;

for func_num = 1:15

    %for overlapping function D = 905
    if func_num > 12 && func_num < 15
        D = 905; % dimensionality of the objective function.
    else
        D = 1000;
    end

    % Search Range
    if (func_num == 1 | func_num == 4 | func_num == 7 | func_num == 8 | func_num == 11 | func_num == 12 | func_num == 13 | func_num == 14 | func_num == 15 )
        lb = -100;
        ub = 100;
    end
    if (func_num == 2 | func_num == 5 | func_num == 9)
        lb = -5;
        ub = 5;
    end
    if (func_num == 3 | func_num == 6 | func_num == 10)
        lb = -32;
        ub = 32;
    end

   for run = 1:runs
      initial_flag = 0; % should set the flag to 0 for each run, each function
      pop = lb + (ub-lb)*rand(D, NP);
      val = benchmark_func(pop, func_num); % fitness evaluation
      FEs = NP;
      while (FEs <= Max_FEs-NP)
         % random search, you should add your own method
         newpop = lb + (ub-lb)*rand(D, NP);
         newval = benchmark_func(newpop, func_num);
         FEs = FEs + NP;
         
         index = (newval < val);
         val(index) = newval(index);
         pop(:, index) = newpop(:, index);
         
         % demo output, you should save your results to some files
         fprintf(1, 'func_num = %d, run = %d, FEs = %d\n', func_num, run, FEs);
         fprintf(1, 'min(val) = %e\n\n', min(val));
      end
   end

end

toc;
