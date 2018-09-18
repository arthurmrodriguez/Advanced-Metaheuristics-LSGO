function [PopOld,Fit,Par]= CC_ANDE(CC_nfes,Pop,Fit,lu,func_num,CC_Alg_Ind,Par)
global initial_flag
initial_flag=initial_flag;

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
    
    R = Gen_R(NP,3);
    R(:,1)=[];
    fr=Fit(R);
    [B,I] = sort(fr,2);
    R_S=[];
    for i=1:NP
        R_S(i,:)=R(i,I(i,:));
    end
    rb=R_S(:,1);
    rm=R_S(:,2);
    rw=R_S(:,3);
    
    F = 0.20+0.6*rand(NP,D);
    %     F = repmat(F,1,D);
    
    p1 = ones(NP,1);
    p2 = 0.75 + 0.25*rand(NP,1);
    p3 = 0.50 + 0.25*rand(NP,1);
    
    p=[p1 p2 p3];
    w=p./repmat(sum(p,2),1,3);
    
    w1= repmat(w(:,1),1,D);
    w2= repmat(w(:,2),1,D);
    w3= repmat(w(:,3),1,D);
    
    
    X = w1.*Pop(rb, :)+w2.*Pop(rm, :)+ w3.*Pop(rw, :)...
        + 2*F.*(Pop(rb, :) - Pop(rw, :));
    
    
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
    
%     fprintf('FES_Perc:%.2f%% FES_No:%1.3e\tFitiBest:%1.3e Tri\n',(Cnfes+nfes)*100/3.00E+06,Cnfes+nfes,min(Fit));
    
end

PopOld(:,CC_Alg_Ind)=Pop;
Par.stgCount=stgCount;
Par.CRNewFlags=CRNewFlags;
Par.CR=CR;
Par.CRRatio=CRRatio;

end