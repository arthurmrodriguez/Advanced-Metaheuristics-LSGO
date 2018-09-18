function [A,CR,stgCount]=Cr_Adaptation(CRNewFlags,GenRatio,stgCount,CRRatio,CRs)

if(GenRatio<=(1/10))
    if(GenRatio<=(1/60))
        A=[0.05 0.1];
    elseif(GenRatio<=(1/40) && GenRatio>(1/60))
        A=[0.05 0.1 0.2 0.3];
    elseif(GenRatio<=(1/30) && GenRatio>(1/40))
        A=[0.05 0.1 0.2 0.3 0.4 0.5];
    elseif(GenRatio<=(1/24) && GenRatio>(1/30))
        A=[0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7];
    elseif(GenRatio<=(1/20) && GenRatio>(1/24))
        A=[0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
    elseif(GenRatio<=(1/10) && GenRatio>(1/20))
        A=[0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
    end
else
    A=[0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.95];
end


CR=CRs;
CR_New_ind=(CRNewFlags==0);

if(sum(CR_New_ind)>0)
    if(GenRatio<=(1/10))
        paraIndex=ceil(length(A)* rand(sum(CR_New_ind), 1));
        CR(CR_New_ind)=A(paraIndex);
    else
        stgCount(CR_New_ind)=stgCount(CR_New_ind)+1;
        stgCount_ind=(stgCount==21);
        paraIndex=ceil(length(A)* rand(sum(stgCount_ind), 1));
        stgCount(stgCount_ind)=0;
        CR(stgCount_ind)=A(paraIndex);
    end
end

CR_Imp_ind=(CRNewFlags==1);

if(sum(CR_Imp_ind)>0)
    CR(CR_Imp_ind)=A(abs(CRRatio)==max(abs(CRRatio)));
end


end

