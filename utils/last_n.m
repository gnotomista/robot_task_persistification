function v_n = last_n(v, n)
if numel(v) > n
    v_n = v(end-n+1 : end);
else
    v_n = v;
end
end