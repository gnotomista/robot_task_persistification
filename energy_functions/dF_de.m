function dF_de_out = dF_de(x, e, t, k, lambda, Ic)

dF_de_out = k * dw_de(x, e, t, lambda, Ic) - k;

end