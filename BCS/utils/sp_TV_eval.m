% evaluates total variation (TV) norm  

function [out]= sp_TV_eval(z)
    tv_sp = sqrt(abs(z{1}).^2 + abs(z{2}).^2);
    out = sum(abs(tv_sp(:)));
  end