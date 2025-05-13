function[P,indices] = calc_deim_proj_mat(U_nl)
    [maxval, ind] = max(abs(U_nl(:,1)), [], 1);
    [n,m] = size(U_nl);
    P = zeros(n,m);
    indices = zeros(1,m);
    P(ind,1) = 1;
    indices(1,1) = ind;
    for i=2:m;
    % These two ways of calculating c should give the same result
        %c = (P(:,1:i-1)'*U_nl(:,1:i-1)) \ P(:,1:i-1)'*U_nl(:,i)
        c = U_nl(indices(1:i-1),1:i-1) \ U_nl(indices(1:i-1),i);
    r = U_nl(:,i) - U_nl(:,1:i-1)*c;
    [maxval, ind] = max(abs(r), [], 1);
    indices(1,i) = ind;
    P(ind,i) = 1;
    end;
end

