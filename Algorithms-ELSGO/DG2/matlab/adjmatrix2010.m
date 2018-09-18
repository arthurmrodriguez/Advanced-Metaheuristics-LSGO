
% Number of non-separable groups for each function in CEC'210 benchmark suite.
numNonSep = [0 0 0 1 1 1 1 1 10 10 10 10 10 20 20 20 20 20 1 1];
 
for f=[1:1:20]
    f
    filename1 = sprintf('./cec2010/datafiles/f%02d_op.mat', f);
    filename2 = sprintf('./cec2010/datafiles/f%02d_opm.mat', f);
    flag = false;
    if(exist(filename1))
        load(filename1);
        flag = true;
    elseif(exist(filename2))
        load(filename2);
        flag = true;
    end

    dim = 1000;
    Adj = zeros(dim,dim);

    m = 50;
    if(flag)
        for g=1:numNonSep(f)
            lb = (g-1)*m+1;
            ub = g*m;
            l = p((g-1)*m+1:g*m);
            [a b] = meshgrid(l, l);
            pairs = [a(:) b(:)];
            size(pairs)

            if(ismember(f, [8, 13, 18]))
                for q=1:length(l)-1
                    c1 = min([l(q), l(q+1)]);
                    c2 = max([l(q), l(q+1)]);
                    Adj(c1, c2) = 1;
                    Adj(c2, c1) = 1;
                end
            else
                for k=1:size(pairs, 1)
                    coord = pairs(k, :);
                    Adj(coord(1), coord(2)) = 1;
                    Adj(coord(2), coord(1)) = 1;
                end
            end


        end
    end

    if(f == 19)
        Adj = ones(1000);
    end
    if(f == 20)
        Adj = diag(ones(1,999), 1);
        Adj = Adj + Adj';
    end
    Adj = Adj | diag(diag(eye(dim)));
    filename = sprintf('./cec2010/adjacency/f%02d.mat', f);
    save('-7', filename, 'Adj');

end

