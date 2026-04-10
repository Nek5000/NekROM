function[P,indices] = calc_qdeim_proj_mat(U_nl)
    [A,B,P] = qr(U_nl');
    P = P(:,1:size(U_nl,2));
    indices = [];
    for col=1:size(P,2);
        [val, ind] = max(P(:,col));
        indices = [indices, ind];
    end;
end
