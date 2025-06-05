function[c] = lcurl(u,v,x,y)
    uy=lgrad(u,x,y,2);
    vx=lgrad(v,x,y,1);
    c = vx-uy;
end

