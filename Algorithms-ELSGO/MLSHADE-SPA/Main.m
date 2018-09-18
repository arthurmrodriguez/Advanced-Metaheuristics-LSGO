%%%%%%%%%%%%%%%%%%%
%% This package is a MATLAB/Octave source code of MLSHADE-SPA which is an improved version of LSHADE-SPA.
%% Please see the following paper:
%% (Under publication) Anas A. Hadi, Ali W. Mohamed, and Kamal M. Jambi: LSHADE-SPA Memeteic Framework for Solving Large Scale Problems

%% About LSHADE-SPA, please see following papers:
%% Ali W. Mohamed, Anas A. Hadi, Anas M. Fattouh, and Kamal M. Jambi: L-SHADE with Semi Parameter Adaptation Approach for Solving CEC 2017 Benchmark Problems, Proc. IEEE Congress on Evolutionary Computation (CEC-2017), Spain, June, 2017
%% Ryoji Tanabe and Alex Fukunaga: Improving the Search Performance of SHADE Using Linear Population Size Reduction,  Proc. IEEE Congress on Evolutionary Computation (CEC-2014), Beijing, July, 2014.
%% J. Zhang, A.C. Sanderson: JADE: Adaptive differential evolution with optional external archive,” IEEE Trans Evol Comput, vol. 13, no. 5, pp. 945–958, 2009

clear all
clc
%max_nfes = 3.00E+06; % number of generations 5000 for 1.25e+05 FES, 12000 FOR 6.00E+05 FES AND 60000 FOR 3.00E+06 FES
max_nfes = 5.00E+04;
runs=1;%25;
Alg_Name='MLSHADE-SPA';

% for func_num=1:15
%     if func_num > 12 && func_num < 15
%         D = 905;
%     else
%         D = 1000;
%     end
%     if (func_num == 1 || func_num == 4 || func_num == 7 || func_num == 8 || func_num == 11 || func_num == 12 || func_num == 13 || func_num == 14 || func_num == 15 )
%         lb = -100;
%         ub = 100;
%     end
%     if (func_num == 2 || func_num == 5 || func_num == 9)
%         lb = -5;
%         ub = 5;
%     end
%     if (func_num == 3 || func_num == 6 || func_num == 10)
%         lb = -32;
%         ub = 32;
%     end

    % Start time
    time_start = cputime;
    
    % Negative func_num indicates noise levels
    func_num = -256;
    D = 1024;
    dim = abs(D/func_num);
    lb = -8;
    ub = 8;
    problem_id='D';
    % Filename to save data from experiments
    problem_id = strcat(problem_id, int2str(dim));
    if(func_num < 0)
        problem_id = strcat(problem_id,'N');
    end
    output_file = strcat('Results/MLSHADE-SPA_', problem_id,'_', int2str(max_nfes), '.txt');
    
    % Get the bsf_solution to display at the end
    bsf_solution = 0;
    outcome=[];
    Conv_Fit=[];
    fprintf('\n-------------------------------------------------------\n')
    fprintf('Running %s on Function = %s, Dimension size = %d\n',Alg_Name, problem_id, D)
    for run=1:runs
        
            % Prepare arguments
            [bsf_solution, f, All_Fit] = MLSHADE_SPA(max_nfes,lb,ub,func_num,D);
            outcome = [outcome f];
            fprintf('run:%d \t\t FitiBest: %e\n',run ,f);
            file_name=sprintf('Figures\\%s_CEC2013_Problem#%s_Run#%s',Alg_Name,int2str(func_num),int2str(run));
            save(file_name,'All_Fit');
    end
    
    % Get elapsed time
    time_end = cputime-time_start;
    
    % Write results in extern file
    fileID =  fopen(output_file,'w');
    fprintf(fileID,'FitiBest: %e\n',f);
    fprintf(fileID,'Min: %e\n',min(outcome));
    fprintf(fileID,'Max: %e\n',max(outcome));
    fprintf(fileID,'Median: %e\n',median(outcome));
    fprintf(fileID,'Mean: %e\n',mean(outcome));
    fprintf(fileID,'Std: %e\n',std(outcome));
    fprintf(fileID, 'Time: %f\n',time_end);
    fprintf(fileID,'\nBEST:\n');
    fprintf(fileID,'%f ', bsf_solution);
    fprintf(fileID,'\n');
    fclose(fileID);
    
    file_name=sprintf('Results\\%s_CEC2013_Problem#%s',Alg_Name,int2str(func_num));
    save(file_name,'outcome');
%end