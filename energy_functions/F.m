function F_out = F(x, e, t, k, lambda, Ic)

F_out = k * (w(x, e, t, lambda, Ic) - e);

end