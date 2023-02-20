function dF_dx_out = dF_dx(x, e, t, k, lambda, Ic)

dF_dx_out = k * dw_dx(x, e, t, lambda, Ic);

end