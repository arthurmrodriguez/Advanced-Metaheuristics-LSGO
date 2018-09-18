
numNonSep = [0 0 0 7 7 7 7 20 20 20 20 1 20 20 1] ;

for i=[13:1:15]
    filename = sprintf('./cec2013/datafiles/f%02d.mat', i);
    i
    m = 0;
    dim = 1000;
    if(ismember(i, [13 14]))
        m = 5;
        dim = 905;
    end
    Adj = zeros(dim,dim);
    p = 1:dim;
    s = [dim];
    if(exist(filename))
        load(filename);
    end

    c = cumsum(s);

    for g=1:length(s)
        if g == 1
            ldim = 1;
        else
            ldim = c(g-1) - ((g-1)*m) + 1;
        end 
        udim = c(g) - ((g-1)*m);

        l = p(ldim:udim);
        [a b] = meshgrid(l, l);
        pairs = [a(:) b(:)];

        for k=1:size(pairs, 1)
            coord = pairs(k, :);
            Adj(coord(1), coord(2)) = 1;
            Adj(coord(2), coord(1)) = 1;
        end
    end

    if(i == 12)
        Adj = diag(ones(1,999), 1);
        Adj = Adj + Adj';
    end
    Adj = Adj | diag(diag(eye(dim)));
    filename = sprintf('./cec2013/adjacency/f%02d.mat', i);
    save('-7', filename, 'Adj');
end
