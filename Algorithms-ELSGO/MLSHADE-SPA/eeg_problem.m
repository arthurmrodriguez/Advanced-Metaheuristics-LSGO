% EEG Problem objective function implementation
% Adaptation for MLSHADE-SPA algorithm

% EEG problem function: evaluates each genome of the population
% and returns a matrix of evaluations
function Fit = eeg_problem(population,noise)

    aux_pop = population';
    
    Fit = zeros(1,size(aux_pop,1));
    for i = 1:size(aux_pop,1)
        Fit(i) = objective(aux_pop(i,:),noise);
    end
    
end


% Objective function to be evaluated
function fit = objective(genome, noise)
    
    persistent problem_name;
    persistent Mixed ICA_component A_matrix optimum_s1;
    % Update problem_name according to dim
    dim = size(genome,2);
    
    % Initialize problem_name for further evals
    if(isempty(problem_name))
        if(dim == 1024)
            problem_name = "D4";
        elseif(dim == 3072)
            problem_name = "D12";
        else
            problem_name = "D19";
        end
        
        % Check if noise
        if(noise)
           problem_name = problem_name + "N";
        end
        
        % DATAin global Path
        path_to_datain = "../DATAin/";
        
        % Load each matrix only ONCE after problem_name
        % has been established. Persistent
        % values are retained in memory between 
        % calls to the function
        A_matrix = importdata(path_to_datain+problem_name+"A.txt",' ');
        ICA_component = importdata(path_to_datain+problem_name+"S.txt",' ');
        Mixed = importdata(path_to_datain+problem_name+"X.txt",' ');    
        optimum_s1 = importdata(path_to_datain+problem_name+"S1.txt", ' ');
    end
    
    % Objective function
    % Turn array into matrix: reshape appends elements column-wise
    % By shifting dimensions (cols, rows) we ensure S1 has elements in
    % order. Then transpose matrix is used to reestablish original dim
    matrix_S1 = reshape(genome, [size(ICA_component,2),size(ICA_component,1)]);
    matrix_S1 = matrix_S1';
    
    % Calculate X1 = AxS1
    X1 = A_matrix*matrix_S1;
    
    % Calculate COR1 = covar(X,X1)
    COR1 = my_correlation(X1,Mixed);
    
    % Calculate f2
    musum = (ICA_component - matrix_S1).^2;
    musum = sum(sum(musum));
    
    % Fitness value for Single-Objective problem: MIN(f1+f2)
    fit = my_diagonal(COR1) + musum/(size(ICA_component,1)*size(ICA_component,2));
end


% Function to obtain the mean and stdev of 
% a set of values
function [v_mean,v_stdev] = new_mean_std(v)
    v_mean = mean(v);
    v_stdev = std(v);
end

% Function to calculate the COVAR of elements
function COVAR = my_correlation2(A, B)

    COVAR = 0;
    % Mean & stdev of each array
    [A1_mean,A1_stdev] = new_mean_std(A);
    [B1_mean,B1_stdev] = new_mean_std(B);

    a = A1_stdev*B1_stdev;
    
    if(abs(a) > 0.00001)
        temp1 = A - A1_mean;
        temp2 = B - B1_mean;
        COVAR = sum(temp1.*temp2)/(size(A,2)*a);
    end
   
end

% Function to create the effective Covariance matrix of
% between matrix X and AxS1
function COVAR_M = my_correlation(A, B)
    % Preallocate memory for speed
    COVAR_M = zeros(size(A,1), size(A,1));
    for i= 1:size(A,1) 
        for j = 1:size(B,1)
            COVAR_M(i,j) = my_correlation2(A(i,:),B(j,:));
        end
    end
end


% Function to obtain f1, where 
% 1. Non-diagonal elements: 1/(N^2 - N) SUM(i) SUM(j!=i) C_ij^2
% 2. Diagonal elements: 1/N SUM(i) 1-(C_ii^2)

function f1_value = my_diagonal(M)

    % Non-diagonal elements will be squared
    % First we substract them from M and then operate
    non_diagonal = M - diag(diag(M));
    partial_non_diagonal = sum(sum(non_diagonal.^2)) / (size(M,1)^2 - size(M,1));
    
    % Diagonal elements will be operated differently
    partial_diagonal = sum(sum((1-(diag(M))).^2)) / size(M,1);
    
    f1_value = partial_non_diagonal + partial_diagonal;
end    