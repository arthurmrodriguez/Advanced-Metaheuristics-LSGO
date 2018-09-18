function r=genR_EADE(Fit)

NP=length(Fit);
r(:,1)=1:NP;

[srt, Fit_index]=sort(Fit,'ascend');

T=ceil(length(Fit_index)/10);
Best=Fit_index(1:T);
Mid=Fit_index(T+1:end-T);
Worest=Fit_index(end-T+1:end);

% choose three random individuals from population mutually different
r(:,2) = Best(ceil(length(Best)* rand(NP, 1)));

r(:,3) = Worest(ceil(length(Worest)* rand(NP, 1)));

r(:,4) = Mid(ceil(length(Mid)* rand(NP, 1)));

pos=r(:,2)==r(:,3);

while(sum(pos)~=0)
    r(pos,3) = Worest(ceil(length(Worest)* rand(sum(pos), 1)));
    pos=r(:,2)==r(:,3);
end

pos=r(:,3)==r(:,4);
while(sum(pos)~=0)
    r(pos,4) = Mid(ceil(length(Mid)* rand(sum(pos), 1)));
    pos=r(:,3)==r(:,4);
end

end

