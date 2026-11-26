function val = get_opt(opts, name, defaultVal)
%GET_OPT  Field of an options struct, or a default value if it is missing.
%
%   val = GET_OPT(opts, name, defaultVal)
%
%   Example:
%     opts.maxIter = 20;
%     a = get_opt(opts, 'tol', 1e-6);      % a = 1e-6
%     b = get_opt(opts, 'maxIter', 100);   % b = 20
    if isfield(opts, name)
        val = opts.(name);
    else
        val = defaultVal;
    end
end
