function [X_best,Best_fit,All_fit]= MLSHADE_SPA(max_nfes,L,H,func_num,D)

global initial_flag
initial_flag=0;
format long

lu = [L * ones(1, D); H * ones(1, D)];

NP=250;

max_NP = NP;
min_NP = 20.0;

Pop = repmat(lu(1, :), NP, 1) + rand(NP, D) .* (repmat(lu(2, :) - lu(1, :), NP, 1));

Fit = benchmark_func(Pop',func_num);

Fit=Fit';
[~,I_best]=min(Fit);
X_best=Pop(I_best,:);

CC_nfes=max_nfes/50;

Par=[];

Par.stgCount=zeros(NP,1);
Par.CRNewFlags=zeros(NP,1);
Par.CR=zeros(NP,1);
Par.CRRatio=zeros(11,1);

Par.memory_size=5;
Par.memory_sf = 0.5 .* ones(Par.memory_size, 1);
Par.memory_pos = 1;

archive=[];
archive.NP = NP; % the maximum size of the archive
archive.Pop = []; % the solutions stored in te archive;
archive.funvalues = []; % the function value of the archived solutions

nfes = NP;

flag=0;

All_fit=[nfes,min(Fit)];

while nfes<max_nfes
    
    Par.GenRatio=nfes/max_nfes;
    
    EA_nfes =round(CC_nfes/2);
    EA_nfes=EA_nfes-mod(EA_nfes,NP);
    
    MMTS_nfes =round(CC_nfes/2);
    MMTS_nfes=MMTS_nfes-mod(MMTS_nfes,NP);
    
%% EA part
    CC_Group_Ind=ones(1,D);
    % make allocated nfes dividable by NP
    Alg_fit=round(0.5*EA_nfes)-mod(round(0.5*EA_nfes),NP);
    if Alg_fit+nfes>max_nfes
        Alg_fit=max_nfes-nfes;
    end
    [Pop,Fit,archive,Par]= CC_LSHADESPA(Alg_fit,Pop,Fit,lu,func_num,CC_Group_Ind==1,archive,Par);
    
    nfes = nfes+Alg_fit;
    
    if nfes >= max_nfes;
        [Best_fit,I_best]= min(Fit);
        X_best=Pop(I_best,:);
        All_fit=[All_fit;nfes,Best_fit];
        return;
    end
    
    All_fit=[All_fit;nfes,min(Fit)];

    
    if(flag==0) % First round
        Alg_CC_nfes=(0.5*EA_nfes/3);
        LSHADESPA_CC_nfes =round(Alg_CC_nfes);
        ANDE_CC_nfes =round(Alg_CC_nfes);
        EADE_CC_nfes =0.5*EA_nfes-(LSHADESPA_CC_nfes+ANDE_CC_nfes);
        flag=1;
    else
        % 90% for previous value and 10% for new one
        LSHADESPA_CC_nfes =round(0.9*LSHADESPA_CC_nfes+0.1*0.5*EA_nfes*All_Imp(1));
        ANDE_CC_nfes =round(0.9*ANDE_CC_nfes+0.1*0.5*EA_nfes*All_Imp(2));
        EADE_CC_nfes =0.5*EA_nfes-(LSHADESPA_CC_nfes+ANDE_CC_nfes);
    end
    
    Group_No=3;
    CC_Group_Ind=ceil(Group_No*rand(1,D));
    while(length(unique(CC_Group_Ind))~=Group_No)
        CC_Group_Ind=ceil(Group_No*rand(1,D));
    end
    
    Alg_Group_Ind = (CC_Group_Ind==1);
    % make allocated nfes dividable by NP
    LSHADESPA_CC_nfes=LSHADESPA_CC_nfes-mod(LSHADESPA_CC_nfes,NP);
    if (LSHADESPA_CC_nfes+nfes)>max_nfes
        LSHADESPA_CC_nfes=max_nfes-nfes;
    end
    [Pop,LSHADESPA_Fit,archive,Par]= CC_LSHADESPA(LSHADESPA_CC_nfes,Pop,Fit,lu,func_num,Alg_Group_Ind,archive,Par);
    nfes = nfes+LSHADESPA_CC_nfes;
    
    if nfes >= max_nfes;
        [Best_fit,I_best]= min(Fit);
        X_best=Pop(I_best,:);
        All_fit=[All_fit;nfes,Best_fit];
        return;
    end
    
    Alg_Group_Ind = (CC_Group_Ind==2);
    % make allocated nfes dividable by NP
    ANDE_CC_nfes=ANDE_CC_nfes-mod(ANDE_CC_nfes,NP);
    if (ANDE_CC_nfes+nfes)>max_nfes
        ANDE_CC_nfes=max_nfes-nfes;
    end
    [Pop,ANDE_Fit,Par]= CC_ANDE(ANDE_CC_nfes,Pop,LSHADESPA_Fit,lu,func_num,Alg_Group_Ind,Par);
    nfes = nfes+ANDE_CC_nfes;
    
    if nfes >= max_nfes;
        [Best_fit,I_best]= min(Fit);
        X_best=Pop(I_best,:);
        All_fit=[All_fit;nfes,Best_fit];
        return;
    end
    
    Alg_Group_Ind = (CC_Group_Ind==3);
    % make allocated nfes dividable by NP
    EADE_CC_nfes=EADE_CC_nfes-mod(EADE_CC_nfes,NP);
    if (EADE_CC_nfes+nfes)>max_nfes
        EADE_CC_nfes=max_nfes-nfes;
    end
    [Pop,EADE_Fit,Par]= CC_EADE(EADE_CC_nfes,Pop,ANDE_Fit,lu,func_num,Alg_Group_Ind,Par);
    nfes = nfes+EADE_CC_nfes;
    
    if nfes >= max_nfes;
        [Best_fit,I_best]= min(Fit);
        X_best=Pop(I_best,:);
        All_fit=[All_fit;nfes,Best_fit];
        return;
    end
    
    All_fit=[All_fit;nfes,min(Fit)];
    
    %% Calculate Improvmant ratio for each EA
    All_Imp=[];
    All_Imp(:,1)=Fit-LSHADESPA_Fit;
    All_Imp(:,2)=LSHADESPA_Fit-ANDE_Fit;
    All_Imp(:,3)=ANDE_Fit-EADE_Fit;
    All_Imp=sum(All_Imp);
    All_Imp=All_Imp./NP;
    if(max(All_Imp)~=0)
        All_Imp=All_Imp./[LSHADESPA_CC_nfes,ANDE_CC_nfes,EADE_CC_nfes];
        All_Imp=All_Imp./sum(All_Imp);
        [temp_imp,Imp_Ind] = sort(All_Imp);
        for imp_i=1:length(All_Imp)-1
            All_Imp(Imp_Ind(imp_i))=max(All_Imp(Imp_Ind(imp_i)),0.1); % No one work lessthan 10%
        end
        All_Imp(Imp_Ind(end))=1-sum(All_Imp(Imp_Ind(1:end-1)));
    else
        Imp_Ind=1:length(All_Imp);
        All_Imp(:)=1/length(All_Imp);
    end
    Fit=EADE_Fit;

    %% MMTS
    CC_Group_Ind=ones(1,D);
    MMTS_Group_Ind = (CC_Group_Ind==1);
    if MMTS_nfes+nfes>max_nfes
        MMTS_nfes=max_nfes-nfes;
    end
    [Pop,Fit,Par]= MMTS(MMTS_nfes,Pop,Fit,lu,func_num,MMTS_Group_Ind,Par);
    
    nfes=nfes+MMTS_nfes;
    if nfes >= max_nfes;
        [Best_fit,I_best]= min(Fit);
        X_best=Pop(I_best,:);
        All_fit=[All_fit;nfes,Best_fit];
        return;
    end
    
    %% Population size reduction 
    plan_NP = round((((min_NP - max_NP) / (0.5*max_nfes)) * nfes) + max_NP); % reach min in 1/2 max_nfes
    
    if NP > plan_NP
        reduction_ind_num = NP - plan_NP;
        if NP - reduction_ind_num <  min_NP; reduction_ind_num = NP - min_NP;end
        NP = NP - reduction_ind_num;
        for r = 1 : reduction_ind_num
            [valBest, indBest] = sort(Fit, 'ascend');
            worst_ind = indBest(end);
            Pop(worst_ind,:) = [];
            Fit(worst_ind,:) = [];
            Par.stgCount(worst_ind)=[];
            Par.CRNewFlags(worst_ind)=[];
            Par.CR(worst_ind)=[];
        end
        archive.NP = NP;
        if size(archive.Pop, 1) > archive.NP
            rndpos = randperm(size(archive.Pop, 1));
            rndpos = rndpos(1 : archive.NP);
            archive.Pop = archive.Pop(rndpos, :);
        end
    end
    
    fprintf('FES_Perc:%.2f%%\tFES:%1.3e\tFitiBest:%1.3e\n',nfes*100/max_nfes,nfes,min(Fit));
    
    
end

[Best_fit,I_best]= min(Fit);
X_best=Pop(I_best,:);
All_fit=[All_fit;nfes,Best_fit];