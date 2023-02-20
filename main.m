clc
clear
close all

% constants
DT = 0.02;
T_MAX = 1e3;
X_MIN = -2;
X_MAX = 2;
% U_MAX = 1;
OPTIM_OPTIONS = optimoptions(@quadprog, 'Display', 'off');
N_TRAJ_PLOT = 200;
X_VEC= linspace(X_MIN, X_MAX, 40);
Y_VEC = linspace(X_MIN, X_MAX, 40);
[XX, YY] = meshgrid(X_VEC, Y_VEC);
ZZ = zeros(size(XX));
N_LEVELS = 10;
% OPTION 1: $h_1 = -\dot V \ge 0$
% OPTION 2: $h_1 = -\dot V - gamma V_1 \ge 0$
OPTION = 2;

% task definition (u_hat)
p = @(t) 1 * [cos(t); sin(t)];
pdot = @(t) 1 * [-sin(t); cos(t)];
u_hat = @(x, t) pdot(t) + p(t) - x;

% energy
e_chg = 0.9;
e_min = 0.6;
k = 3;
lambda = 0.015;
Ic = 0.85;

% cbf-clf-qp
switch OPTION
    case 1
        gamma1 = 1e0;
        gamma2 = 1e-1;
        gamma_V1 = 1e-2;
    case 2
        gamma1 = 1e0;
        gamma2 = 1e-1;
        gamma_V1 = 1e-2;
        gamma_V2 = 1e-1;
end

% state initialization
x = [1; 1]; % X_MIN + (X_MAX - X_MIN) * rand(2, 1);
e = 0.7;

% plot initialization
figure('units', 'normalized', 'position', [0 0 1 1])
% energy
subplot(3, 2, 1), hold on, grid on
ylim([0 1])
he = plot(0, e, 'Color', [0.1 0.8 0.5], 'LineWidth', 2);
title('energy')
% input difference
subplot(3, 2, 2), hold on, grid on
hc = plot(0, nan, 'Color', [0.7 0.3 0.1], 'LineWidth', 2);
title('input difference')
% workspace
subplot(3, 2, 3 : 6), hold on, grid on, axis equal, axis([X_MIN X_MAX X_MIN X_MAX])
hC = plot_environment(0, Ic, [], X_VEC, Y_VEC, XX, YY, ZZ, N_LEVELS);
hp = scatter([1 0] * p(0), [0 1] * p(0), 500, [0.8 0.2 0.3], '.');
hrt = plot(x(1), x(2), 'Color', [0.2 0.4 0.6], 'LineWidth', 2);
lgd = legend([hp, hrt, hC{1}],'reference', 'robot trajectory', 'intensity level curves ($I = I_c$ in green)');
lgd.Interpreter = 'latex';
lgd.FontSize = 12;
title('workspace')

for t = 0 : DT : T_MAX
    tic

    F_value = F(x, e, t, k, lambda, Ic);
    dF_dx_value = dF_dx(x, e, t, k, lambda, Ic);
    dF_de_value = dF_de(x, e, t, k, lambda, Ic);
    dF_dt_value = dF_dt(x, e, t, k, lambda, Ic);

    e_dot = F_value;

    % cbf
    h1 = (e_chg - e) * (e - e_min);
    h1_dot = (e_chg + e_min - 2 * e) * F_value;
    h2 = h1_dot + gamma1 * h1;
    Acbf = - (e_chg + e_min - 2 * e) * dF_dx_value;
    bcbf = (e_chg + e_min - 2 * e) * dF_dt_value + ...
        (-2 * F_value + (e_chg + e_min - 2 * e) * dF_de_value + ...
        gamma1 * (e_chg + e_min - 2 * e)) * F_value + ...
        gamma2 * h2;

    % clf
    switch OPTION
        case 1
            V = (1 - e) ^ 2;
            V_dot = -2 * (1 - e) * F_value;
            h1 = -V_dot;
            Aclf = -2 * (1 - e) * dF_dx_value;
            bclf = 2 * (1 - e) * dF_dt_value + ...
                (-2 * F_value + 2 * (1 - e) * dF_de_value) * F_value + ...
                gamma_V1 * h1;
        case 2
            V = (1 - e) ^ 2;
            V_dot = -2 * (1 - e) * F_value;
            h1 = -V_dot - gamma_V1 * V;
            Aclf = -2 * (1 - e) * dF_dx_value;
            bclf = (-2 * F_value + 2 * (1 - e) * dF_de_value - 2 * gamma_V1 * (1 - e)) * F_value + ...
                2 * (1 - e) * dF_dt_value - 2 * gamma_V1 * (1 - e) * F_value + ...
                gamma_V2 * h1;
    end

    % QP
    u_hat_value = u_hat(x, t);
    kappa_value = kappa(e, e_dot, e_min, e_chg);
    Hqp = blkdiag(2 * eye(2), 2 * kappa_value);
    fqp = [-2 * u_hat_value; 0];
    Aqp = [Acbf 0; Aclf -1];
    bqp = [bcbf; bclf];
    u_delta = quadprog(Hqp, fqp, Aqp, bqp, [], [], [], [], [], OPTIM_OPTIONS);

    % integration step
    u = u_delta(1:2);
    % if norm(u) >= U_MAX
    %     u = u / norm(u) * U_MAX;
    % end
    x = x + u * DT;
    e = e + e_dot * DT;

    % print
    % disp({'e, e_dot, kappa_value', [e, e_dot, kappa_value]})
    % disp({'V_dot', [V_dot]})

    % plot updates
    if mod(t, 1) == 0
        hC = plot_environment(t, Ic, hC, X_VEC, Y_VEC, XX, YY, ZZ, N_LEVELS);
    end
    he.XData = [he.XData, t];
    he.YData = [he.YData, e];
    hc.XData = [hc.XData, t];
    hc.YData = [hc.YData, norm(u-u_hat_value)];
    hrt.XData = last_n([hrt.XData, x(1)], N_TRAJ_PLOT);
    hrt.YData = last_n([hrt.YData, x(2)], N_TRAJ_PLOT);
    hp.XData = [1 0] * p(t);
    hp.YData = [0 1] * p(t);
    drawnow limitrate
    % pause(DT-toc)
end

function hC = plot_environment(t, Ic, hC, X_VEC, Y_VEC, XX, YY, ZZ, N_LEVELS)
for i = 1 : size(XX, 1)
    for j = 1 : size(XX, 2)
        ZZ(i, j) = I_xt([XX(i, j); YY(i, j)], t);
    end
end
xyC = contourc(X_VEC, Y_VEC, ZZ, [Ic, linspace(0.1, 0.99, N_LEVELS)]);
xC = [];
yC = [];
if ~isempty(xyC)
    [xC, yC, zC] = C2xyz(xyC);
end
if isempty(hC)
    hC = cell(numel(xC), 1);
    for i = 1 : numel(xC)
        hC{i} = plot(xC{i}, yC{i}, '-', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2);
    end
else
    for i = 1 : min(numel(xC), numel(hC))
        if zC(i) == Ic
            color = [0.1 0.8 0.1];
            linewidth = 4;
        else
            color = [0.9290 0.6940 0.1250];
            linewidth = 2;
        end
        set(hC{i}, 'XData', xC{i}, 'YData', yC{i}, 'Color', color, 'LineWidth', linewidth)
    end
end
end
