function [phi,phit] = def_phi_phit_ID

  
        phi = @(U) ((U));warning off all;
        phit = @(DU) ((DU));
end