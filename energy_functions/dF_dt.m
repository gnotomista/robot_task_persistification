function dF_dt_out = dF_dt(x, e, t, k, lambda, Ic)

dF_dt_out = k * dw_dt(x, e, t, lambda, Ic);

end