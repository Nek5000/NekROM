function [Me] = get_Me(snaps)
    x=snaps.flds{1}.x;
    y=snaps.flds{1}.y;
    nx1 = size(x,1);
    [zi, w] = zwgll(nx1-1);
    d = deriv_mat(zi);
    [xr,yr,xs,ys,rx,ry,sx,sy,jac,jaci,d] = deriv_geo(x,y,d);
    nL = prod(size(x));
    Me = reshape(jac.*(w*w'),nL,1);
end
