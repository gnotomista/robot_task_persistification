function dw_de_out = dw_de(x, e, t, lambda, Ic)

dw_de_out = 1 / (1 + (1 - e) / e * exp(-lambda * (I_xt(x, t) - Ic)))^2 * ...
    exp(-lambda  * (I_xt(x, t) - Ic)) * 1 / e ^ 2;

end