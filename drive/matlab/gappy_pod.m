% https://epubs.siam.org/doi/pdf/10.1137/22M1484018
function[indices] = gappy_pod(U_nl, n)
    [maxval, ind] = max(abs(U_nl(:,1)), [], 1);
    [N,p] = size(U_nl);
    indices = [ind];
    n_iter = ceil((n-1)/(p-1));
    for i=2:p;
        for k=1:n_iter;
            l=(i-2)*n_iter + k;
            c = pinv(U_nl(indices(1:l),1:i-1), U_nl(indices(1:l),i));
            r = U_nl(:,i) - U_nl(:,1:i-1)*c;
            [maxval, ind] = max(abs(r), [], 1);
            indices(end+1) = ind;
            if len(indices) == n;
                return;
            end;
        end;
    end;
end


