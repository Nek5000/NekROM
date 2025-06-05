function [a, a0, b, c, c0, c1, c2, c3, u0, uk, ukmin, ukmax] =  get_r_dim_ops(au_full, bu_full, cu_full, u0_full, uk_full, nb)
   index  = [1:nb+1];
   index1 = [1:nb];
   index2 = [2:nb+1];
   a      = au_full(index2,index2);
   a0     = au_full(index2,1);
   b      = bu_full(index2,index2);

   cutmp  = cu_full(index1,index,index);
   c0      = reshape(cutmp,nb*(nb+1),nb+1);
   c1 = cutmp(:,1,1);                                                       
   c2 = reshape(cutmp(:,1,:),nb,nb+1);                                      
   c3 = reshape(cutmp(:,:,1),nb,nb+1); 

   c      = cu_full(index1,index2,index2);
   c      = reshape(c,nb*(nb),nb);

   u0      = u0_full(index);

   uk = uk_full(index,:);
   ukmin = min(uk,[],2);
   ukmax = max(uk,[],2);
end

