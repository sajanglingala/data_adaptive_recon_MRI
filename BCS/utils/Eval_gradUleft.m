function left = Eval_gradUleft(U,A,At,phi, phit, V,opts)

U=reshape(U, opts.n1*opts.n2,opts.r);
left = At(A(U*V))*V' + (opts.lambda1*opts.beta1/2)*phit(phi(U));
left = left(:);
