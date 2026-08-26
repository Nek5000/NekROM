function value = read_env_bool(name, default_value)
    % READ_ENV_BOOL Parse a boolean environment variable.
    %
    % Accepts 1/0, true/false, yes/no, on/off. Falls back to default_value
    % when the environment variable is unset.

    env_value = getenv(name);
    if isempty(env_value)
        value = default_value;
        return;
    end

    switch lower(strtrim(env_value))
        case {'1', 'true', 'yes', 'on'}
            value = true;
        case {'0', 'false', 'no', 'off'}
            value = false;
        otherwise
            error('read_env_bool:InvalidValue', ...
                'Invalid boolean value for %s: %s', name, env_value);
    end
end
