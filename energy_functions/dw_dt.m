function dw_dt_out = dw_dt(x, e, t, lambda, Ic)

dw_dt_out = -1 / (1 + (1 - e) / e * exp(-lambda * (I_xt(x, t) - Ic)))^2 * ...
    (-lambda * (1 - e) / e * exp(-lambda * (I_xt(x, t) - Ic))) * dI_dt(x, t);

end