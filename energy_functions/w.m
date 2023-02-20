function w_out = w(x, e, t, lambda, Ic)

w_out = 1 / (1 + (1 - e) / e * exp(-lambda * (I_xt(x, t) - Ic)));

end