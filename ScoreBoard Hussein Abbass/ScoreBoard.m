function Score = ScoreBoard(filename,problemtype,scalefile,MorS)


a = load(filename);


if(MorS=='S')
    if(length(scalefile)>0)
       s = load(scalefile);
       a(:,2) = s(3) + (s(4)-s(3))*(a(:,2)-s(1))/(s(2)-s(1)); 
    end

    baseline = load(strcat('baseline/',problemtype,'S.txt'));
    sol = min(a(:,1));
    if(sol<baseline)
        Score = baseline - sol;
    else
        Score = -1000*(sol - baseline);
    end
    
end
if(MorS=='M')
    if(length(scalefile)>0)
       s = load(scalefile);
       a(:,2) = s(3) + (s(4)-s(3))*(a(:,2)-s(1))/(s(2)-s(1)); 
    end
    
    
    baseline = load(strcat('baseline/',problemtype,'PF.txt'));
    [THETA,RHO] = cart2pol(a(:,1),a(:,2));
    
    W = [a THETA*180/pi];
    
    F1 = sortrows(W,size(W,2));
    
    k=size(F1,1);
    index = [1 round(k/4) round(k/2) round(3*k/4) k];
    F1=F1(index,:);
    Score = 0;
    for i=1:5
	
       if(F1(i,1)>baseline(i,1) && F1(i,2)>baseline(i,2))
           Score = Score - 1000 * ((F1(i,1)-baseline(i,1))^2 + (F1(i,2)-baseline(i,2))^2)^0.5;
       elseif(F1(i,1)<baseline(i,1) && F1(i,2)<baseline(i,2))
           Score = Score + ((F1(i,1)-baseline(i,1))^2 + (F1(i,2)-baseline(i,2))^2)^0.5;
       end
   
    end
    figure
    hold all
    scatter(F1(:,1),F1(:,2),'bo')
    scatter(baseline(:,1),baseline(:,2),'rx')
    legend(sprintf('Your PF: Score %f',Score),'Baseline')
    title ('Pareto Front')
    xlabel('F1')
    ylabel('F2')
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',14)
    
end


