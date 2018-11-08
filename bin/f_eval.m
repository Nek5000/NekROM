function [f] = f_eval(helm,coef,rhs,par,sample_min,sample_max,b)
   bar1 = sum(log(-coef+sample_max)); bar2 = sum(log(coef-sample_min));
   bar = bar1 + bar2;
   f = 0.5*coef'*helm*coef - coef'*rhs - par*( bar) + 0.5*rhs'*inv(helm)*rhs;
end
