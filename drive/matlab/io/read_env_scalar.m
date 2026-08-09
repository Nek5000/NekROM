function value = read_env_scalar(varname, default_value)
% READ_ENV_SCALAR Read a scalar numeric environment variable with fallback.
%
% Returns default_value when the environment variable is unset or cannot be
% parsed into a finite scalar.

    if nargin < 2
        default_value = [];
    end

    raw = getenv(varname);
    if isempty(raw)
        value = default_value;
        return;
    end

    parsed = str2double(strtrim(raw));
    if isempty(parsed) || ~isfinite(parsed)
        value = default_value;
        return;
    end
    value = parsed;
end

