function [PopOld,Fit,Par]= MMTS(CC_nfes,Pop,Fit,lu,func_num,CC_Alg_Ind,Par)
global initial_flag
initial_flag=initial_flag;

Lbound=lu(1,CC_Alg_Ind);
Ubound=lu(2,CC_Alg_Ind);

PopOld=Pop;
Eval_Pop=Pop;
Pop=Pop(:,CC_Alg_Ind);
[NP,D]=size(Pop);

%%%%%% select the MMTS search agents from the DE population by Clearing procedure
[~,in]=sort(Fit);
LS_ind = in(1);
Eval_Pop=Eval_Pop(LS_ind,:);
LS_Pop=Pop(LS_ind,:);
LS_Fit=Fit(LS_ind);
LS_SR=(max(Pop,[],1)-min(Pop,[],1)).*rand(1,D);
LS_SR=min(LS_SR,0.2*(Ubound(1,1:D)-Lbound(1,1:D)));

dim=randperm(D);

nfes=0;

LS_Imp_Flag=1;

while (nfes<=CC_nfes)
    LS_Last_Fit=LS_Fit;
    if LS_Imp_Flag==0
        LS_SR=LS_SR.*rand(1,D);
    end
    for i=1:D
        k=0;
        LS_Flag=1;
        while LS_Flag
            k=k+1;
            LS_Child_pos=LS_Pop;
            LS_Child_pos(dim(i))=LS_Child_pos(dim(i))+k*LS_SR(dim(i));
            if LS_Child_pos(dim(i))>Ubound(i)
                break;
            end
            Eval_Pop(CC_Alg_Ind)=LS_Child_pos;
            LS_Child_fit = benchmark_func(Eval_Pop', func_num);
            nfes=1+nfes;
            if nfes>CC_nfes
                Fit(LS_ind)=LS_Fit;
                Pop(LS_ind,:)=LS_Pop;
                PopOld(:,CC_Alg_Ind)=Pop;
                return;
            end
            if LS_Child_fit<=LS_Fit
                LS_Fit=LS_Child_fit;
                LS_Pop=LS_Child_pos;
            else
                LS_Flag=0;
            end
            
        end
        if k<=1
            k=0;
            LS_Flag=1;
            while LS_Flag;
                k=k+1;
                LS_Child_pos=LS_Pop;
                LS_Child_pos(dim(i))=LS_Child_pos(dim(i))-k*LS_SR(dim(i));
                if LS_Child_pos(dim(i))<Lbound(i)
                    break;
                end
                Eval_Pop(CC_Alg_Ind)=LS_Child_pos;
                LS_Child_fit = benchmark_func(Eval_Pop', func_num);
                nfes=1+nfes;
                if nfes>CC_nfes
                    Fit(LS_ind)=LS_Fit;
                    Pop(LS_ind,:)=LS_Pop;
                    PopOld(:,CC_Alg_Ind)=Pop;
                    return;
                end
                if LS_Child_fit<=LS_Fit
                    LS_Fit=LS_Child_fit;
                    LS_Pop=LS_Child_pos;
                else
                    LS_Flag=0;
                end
                
            end
        end
    end
    if LS_Last_Fit<=LS_Fit
        LS_Imp_Flag=0;
    else
        LS_Imp_Flag=1;
    end
end


Fit(LS_ind)=LS_Fit;
Pop(LS_ind,:)=LS_Pop;
PopOld(:,CC_Alg_Ind)=Pop;
end