function ops_ips = read_ops_ips(ops_dir)
% READ_OPS_IPS Read the offline POD inner-product tag from ops/ips.
%
% Returns H10 if the metadata file is missing or malformed so older cases
% continue to run with the legacy validation default.

    ops_ips = 'H10';
    ips_file = fullfile(ops_dir, 'ips');

    if ~exist(ips_file, 'file')
        warning('read_ops_ips:MissingFile', ...
            'ops/ips not found in %s; defaulting validation to H10.', ops_dir);
        return;
    end

    fid = fopen(ips_file, 'r');
    if fid < 0
        warning('read_ops_ips:OpenFailed', ...
            'Unable to open %s; defaulting validation to H10.', ips_file);
        return;
    end

    raw = fscanf(fid, '%s', 1);
    fclose(fid);

    if isempty(raw)
        warning('read_ops_ips:EmptyFile', ...
            'Empty ops/ips file at %s; defaulting validation to H10.', ips_file);
        return;
    end

    ops_ips = upper(strtrim(raw));
    valid_ops_ips = {'L2', 'H10', 'HLM'};
    if ~ismember(ops_ips, valid_ops_ips)
        warning('read_ops_ips:UnsupportedValue', ...
            'Unsupported ops/ips value "%s"; defaulting validation to H10.', ops_ips);
        ops_ips = 'H10';
    end
end
