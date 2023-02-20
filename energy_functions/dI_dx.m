function dI_dx_out = dI_dx(x, t)

s = 10;
v = 0.25*[1; 1];
sigma = 0.25;
M1 = diag(2 * [-1, sin(0.3 * t / s)]);
M2 = diag(2 * [sin(0.5 * t / s), 1]);

dI_dx_out = exp(-norm(x - M1 * v) ^ 2 / sigma ^ 2) * (-1 / sigma ^ 2) * 2 * (x - M1 * v)' + ...
    exp(-norm(x - M2 * v) ^ 2 / (sigma / 2) ^ 2) * (-1 / (sigma / 2) ^ 2) * 2 * (x - M2 * v)';

end