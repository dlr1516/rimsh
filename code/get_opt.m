% function val = get_opt(opts, name, defaultVal)
%
% Checks if the field 'name' is a member of structure opts.
% If the field is present, then it returns its value.
% Otherwise it returns the default value. 
%
% Example:
%   opts.pippo = 2;
%   a = get_opt(opts, 'pluto', 3.5);  % a = 3.5
%   b = get_opt(opts, 'pippo', -4);   % b = 2
%
function val = get_opt(opts, name, defaultVal)
    if isfield(opts, name)
        val = opts.(name);
    else
        val = defaultVal;
    end
end
