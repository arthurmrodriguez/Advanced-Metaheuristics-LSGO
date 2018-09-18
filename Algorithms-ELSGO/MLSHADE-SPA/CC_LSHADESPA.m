% In this version we test the impact of JADE in our CC framework
function [PopOld,Fit,archive,Par]= CC_LSHADESPA(CC_nfes,Pop,Fit,lu,func_num,CC_Alg_Ind,archive,Par)
global initial_flag
initial_flag=initial_flag;

lu=lu(:,CC_Alg_Ind);

PopOld=Pop;

Eval_Pop=Pop;

Pop=Pop(:,CC_Alg_Ind);

[NP,D]=size(Pop);

X = zeros(NP,D); % trial vector

p_best_rate = 0.1;

stgCount=Par.stgCount;
CRNewFlags=Par.CRNewFlags;
CR=Par.CR;
CRRatio=Par.CRRatio;

memory_size=Par.memory_size;
memory_sf=Par.memory_sf;
memory_pos=Par.memory_pos;

nfes=0;

while nfes<CC_nfes
    
    mem_rand_index = ceil(memory_size * rand(NP, 1));
    mu_sf = memory_sf(mem_rand_index);
    
    [A,CR,stgCount]=Cr_Adaptation(CRNewFlags,Par.GenRatio,stgCount,CRRatio,CR);
    
    [temp_fit, sorted_index] = sort(Fit, 'ascend');
    
    %% for generating scaling factor
    if(nfes <= CC_nfes/2)
        sf=0.45+.1*rand(NP, 1);
        pos = find(sf <= 0);
        
        while ~ isempty(pos)
            sf(pos)=0.45+0.1*rand(length(pos), 1);
            pos = find(sf <= 0);
        end
    else
        sf = mu_sf + 0.1 * tan(pi * (rand(NP, 1) - 0.5));
        
        pos = find(sf <= 0);
        
        while ~ isempty(pos)
            sf(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
            pos = find(sf <= 0);
        end
    end
    sf = min(sf, 1);
    
    r0 = [1 : NP];
    if(size(archive.Pop,1)~=0)
    Arc_pop=archive.Pop(:,CC_Alg_Ind);
    popAll = [Pop; Arc_pop];
    else
        popAll = Pop;
    end
    [r1, r2] = gnR1R2(NP, size(popAll, 1), r0);
    
    pNP = max(round(p_best_rate * NP), 2); %% choose at least two best solutions
    randindex = ceil(rand(1, NP) .* pNP); %% select from [1, 2, 3, ..., pNP]
    randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
    pbest = Pop(sorted_index(randindex), :); %% randomly choose one of the top 100p% solutions
    
    X = Pop+ sf(:, ones(1, D)) .* (pbest - Pop + Pop(r1, :) - popAll(r2, :));

    X = boundConstraint(X, Pop, lu);
    
    mask = rand(NP, D) > CR(:, ones(1, D)); %mask is used to indicate which elements of ui comes from the parent
    Rnd=ceil(D* rand(NP, 1)); %choose one position where the element of X doesn't come from the parent
    jrand = sub2ind([NP D], (1:NP)', Rnd);
    mask(jrand)=false;
    X(mask) = Pop(mask);
                
    Eval_Pop(:,CC_Alg_Ind)=X;
    Child_Fit= benchmark_func(Eval_Pop',func_num);
    Child_Fit=Child_Fit';
    
    nfes = nfes+NP;
    if nfes > CC_nfes;
        PopOld(:,CC_Alg_Ind)=Pop;
        Par.stgCount=stgCount;
        Par.CRNewFlags=CRNewFlags;
        Par.CR=CR;
        Par.CRRatio=CRRatio;
        
        Par.memory_size=memory_size;
        Par.memory_sf=memory_sf;
        Par.memory_pos=memory_pos;
        return;
    end
    
    %% JADE Archive & F Update
    Fit_imp_inf = (Child_Fit<=Fit);

    goodF = sf(Fit_imp_inf);
    dif = abs(Fit - Child_Fit);
    dif_val = dif(Fit_imp_inf);
    
    Eval_Pop(:,CC_Alg_Ind)=Pop;
    archive = updateArchive(archive, Eval_Pop(Fit_imp_inf, :), Fit(Fit_imp_inf));

    num_success_params = numel(goodF);
    if num_success_params > 0
        sum_dif = sum(dif_val);
        dif_val = dif_val / sum_dif;
        
        %% for updating the memory of scaling factor
        memory_sf(memory_pos) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);
        
        memory_pos = memory_pos + 1;
        if memory_pos > memory_size;  memory_pos = 1; end
    end

    %% EADE CR and Pop Update
    CRNewFlags(Fit_imp_inf)=1;
    CRNewFlags(~Fit_imp_inf)=0;
    
    val= Child_Fit./Fit;
    
    val=1-val;
    
    Pop(Fit_imp_inf,:) = X(Fit_imp_inf,:); % replace current by trial
    Fit(Fit_imp_inf) = Child_Fit(Fit_imp_inf) ;
    
    
    for j=1:length(A)
        A_ind=A(j)==CR;
        CRRatio(j)= CRRatio(j)+sum(val(and(A_ind,Fit_imp_inf)));
    end
    
%     fprintf('FES_Perc:%.2f%% FES_No:%1.3e FitiBest:%1.3e SPA\n',(Cnfes+nfes)*100/3.00E+06,Cnfes+nfes,min(Fit));
    
end

PopOld(:,CC_Alg_Ind)=Pop;
Par.stgCount=stgCount;
Par.CRNewFlags=CRNewFlags;
Par.CR=CR;
Par.CRRatio=CRRatio;

Par.memory_size=memory_size;
Par.memory_sf=memory_sf;
Par.memory_pos=memory_pos;
end