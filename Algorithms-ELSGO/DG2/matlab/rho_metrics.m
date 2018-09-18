
function [r1 r2 r3] = rho_metrics(Theta, ThetaI)

    dim = size(ThetaI, 1);

    r1 = sum(sum((Theta .* ThetaI)-eye(dim))) / sum(sum(ThetaI-eye(dim))) * 100;
    r2 = sum(sum(((triu(ones(dim)) - Theta) .* (triu(ones(dim)) - ThetaI)))) / sum(sum(triu(ones(dim)) - ThetaI)) * 100;
    r3 = sum(sum(triu(ones(dim)-abs(Theta-ThetaI))-eye(dim))) * 2 / dim / (dim-1) * 100;

end
