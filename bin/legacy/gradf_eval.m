function gradf = gradf_eval(helm,coef,rhs,par,sample_min,sample_max,b)
   bar = ((-1)./(sample_max-coef) + 1./(coef-sample_min));
   gradf = helm*coef - rhs - par.*bar;
%   gradf = helm'*b*helm*coef - helm*b*rhs - par.*bar;
end
