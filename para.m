function paraconjugate = para(Hz)
%Takes output of dsppset() (symbolic FIR filter) and outputs paraconjugate
    syms z;
    
    paraconjugate = conj(subs(Hz, z, 1/conj(z))); 
end

