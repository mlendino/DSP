function Hz = dsppset(h)
    %This takes in the coefficients for an FIR filter and
    %returns its transfer function
    %I do admit this function is poorly named sorry
    syms z;

    Hz = expand(poly2sym(h,z)/z^(length(h)-1));
    
end

