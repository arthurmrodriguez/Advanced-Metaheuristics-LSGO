% In this version we test the impact of the original EADE as is with in our CC framework
function [PopOld,Fit,Par]= CC_EADE(CC_nfes,Pop,Fit,lu,func_num,CC_Alg_Ind,Par)
global initial_flag
initial_flag=initial_flag;

format long;

lu=lu(:,CC_Alg_Ind);

PopOld=Pop;

Eval_Pop=Pop;

Pop=Pop(:,CC_Alg_Ind);

[NP,D]=size(Pop);


X = zeros(NP,D); % trial vector

stgCount=Par.stgCount;
CRNewFlags=Par.CRNewFlags;
CR=Par.CR;
CRRatio=Par.CRRatio;

nfes=0;

while nfes<CC_nfes
    
    
    [A,CR,stgCount]=Cr_Adaptation(CRNewFlags,Par.GenRatio,stgCount,CRRatio,CR);
    
    mut_prop=rand(NP,1)<=0.5; % Both
    
    r=genR_EADE(Fit);
    
    F1=rand(NP,1);
    F2=rand(NP,1);
    
    X(r(mut_prop,1),:)=Pop(r(mut_prop,4),:) + F1(mut_prop, ones(1, D)).*(Pop(r(mut_prop,2),:)-(Pop(r(mut_prop,4),:))) + F2(mut_prop, ones(1, D)).*((Pop(r(mut_prop,4),:))-(Pop(r(mut_prop,3),:)));
    
    r=Gen_R(NP,4);
    
    temp=sum(~mut_prop);
    
    F=rand(temp,1);
    F(F>0.5)=-1.*F(F>0.5);
    
    F1(~mut_prop)=F;
    
    X(r(~mut_prop,1),:) = Pop(r(~mut_prop,4),:) + F1(~mut_prop, ones(1, D)).* (Pop(r(~mut_prop,2),:) - Pop(r(~mut_prop,3),:));
    
    mask = rand(NP, D) > CR(:, ones(1, D)); %mask is used to indicate which elements of ui comes from the parent
    Rnd=ceil(D* rand(NP, 1)); %choose one position where the element of X doesn't come from the parent
    jrand = sub2ind([NP D], (1:NP)', Rnd);
    mask(jrand)=false;
    X(mask) = Pop(mask);
    
    Temp_Pop = repmat(lu(1, :), NP, 1) + rand(NP, D) .* (repmat(lu(2, :) - lu(1, :), NP, 1));
    
    %% check the lower bound
    xl = repmat(lu(1, :), NP, 1);
    pos = X < xl;
    X(pos)= Temp_Pop(pos);
    
    %% check the upper bound
    xu = repmat(lu(2, :), NP, 1);
    pos = X > xu;
    X(pos)= Temp_Pop(pos);
    
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
        return;
    end
    
    
    Fit_imp_inf = (Child_Fit<=Fit);
    
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
    
    %     fprintf('FES_Perc:%.2f%% FES_No:%1.3e FitiBest:%1.3e EADE\n',(Cnfes+nfes)*100/3.00E+06,Cnfes+nfes,min(Fit));
    
end

PopOld(:,CC_Alg_Ind)=Pop;
Par.stgCount=stgCount;
Par.CRNewFlags=CRNewFlags;
Par.CR=CR;
Par.CRRatio=CRRatio;

end