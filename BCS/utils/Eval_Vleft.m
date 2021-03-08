function out=Eval_Vleft(V,A,At,U,opts)

V=reshape(V,opts.r,opts.n3);
out=2*(U'*At(A(U*V)))+opts.beta2*V;
out=out(:);
end