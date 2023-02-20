function dI_dt_out = dI_dt(x,t)

s = 10;
v = 0.25*[1; 1];
sigma = 0.25;
M1 = diag(2 * [-1, sin(0.3 * t / s)]);
M2 = diag(2 * [sin(0.5 * t / s), 1]);
dM1_dt = diag(2*[0. cos(0.3 * t / s) * 0.3 / s]);
dM2_dt = diag(2*[cos(0.5 * t / s) * 0.5 / s, 0]);

dI_dt_out = exp(-norm(x - M1 * v) ^ 2 / sigma ^ 2) * ...
    (-1 / sigma ^ 2) * 2 * (x - M1 * v)' * (-dM1_dt * v) + ...
    exp(-norm(x - M2 * v) ^ 2 / (sigma / 2) ^ 2) * (-1 / (sigma / 2) ^ 2) * 2 * (x - M2 * v)' * (-dM2_dt * v);

end