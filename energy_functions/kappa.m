function kappa_value = kappa(e, e_dot, e_min, e_chg)

% N = 1e3;
%
% e = (e - e_min) ./ (e_chg - e_min);
% e = min(max(e, 0), 1);
%
% if e_dot > 0
%     kappa_value = nthroot(1 - e .^ N, N);
% else
%     kappa_value = 1 - nthroot(1 - (e - 1) .^ N, N);
% end
%
% kappa_value = 1e6 * kappa_value;

if e_dot > -1e-3
    if e > e_chg - 1e-2
        kappa_value = 0;
    else
        kappa_value = 1e9;
    end
else
    kappa_value = 0;
end

end