function val = pnebi_eval_0(x)

% double t, t2, tinv, val;

if x < 0.0
    x = -x;
end
t = x / 3.75;

if (t < 1.0)
    t2 = t * t;
    val = 1.0 + t2 * (3.5156229 + t2 * (3.0899424 + t2 * (1.2067492 + t2 * (0.2659732 + t2 * (0.360768e-1 + t2 * 0.45813e-2)))));
    val = 2.0 * exp(-x) * val;
else
    tinv = 1 / t;
    val = (0.39894228 + tinv * (0.1328592e-1 + tinv * (0.225319e-2 + tinv * (-0.157565e-2 + tinv * ...
        (0.916281e-2 + tinv * (-0.2057706e-1 + tinv * (0.2635537e-1 + tinv * (-0.1647633e-1 + tinv * 0.392377e-2))))))));
    val = 2.0 * val / sqrt(x);
end

%     return val

end %function
