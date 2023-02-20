function dw_dx_out = dw_dx(x, e, t, lambda, Ic)

dw_dx_out = -1 / (1 + (1 - e) / e * exp(-lambda * (I_xt(x, t) - Ic)))^2 * ...
    (-lambda * (1 - e) / e * exp(-lambda * (I_xt(x, t) - Ic))) * dI_dx(x, t);

end