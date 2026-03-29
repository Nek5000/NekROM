function[ux, uy] = lgrad(u,x,y,mode)
    persistent rx ry sx sy jaci d my_lgrad
    if isempty(rx);
        nx1 = size(x,1);
        [zi, w] = zwgll(nx1-1);
        d = deriv_mat(zi);
        [xr,yr,xs,ys,rx,ry,sx,sy,jac,jaci,d] = deriv_geo(x,y,d);
        my_lgrad=@(u,mode) grad(u,rx,ry,sx,sy,jaci,d,mode);
    end;
    [ux, uy] = my_lgrad(u,mode);
end


