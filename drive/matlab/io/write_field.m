function [] = write_field(basepath, data, wdsize)
    % Determined Dimensions
    sz = size(data.x); 
    ndim = 2; 
    if isfield(data, 'z') || (numel(sz) == 4), ndim = 3; end
    
    nx = sz(1); ny = sz(2);
    nz = 1; if ndim == 3, nz = sz(3); end
    nxyz = nx * ny * nz;

    % 1. Robust Metadata Extraction
    if isfield(data, 'nel'),        nelt = data.nel;
    elseif isfield(data, 'inde'),   nelt = size(data.inde, 1);
    else,                           nelt = sz(end); 
    end

    if isfield(data, 'istep'),      istep_val = data.istep;
    elseif isfield(data, 'iostep'), istep_val = data.iostep;
    else,                           istep_val = 0; 
    end

    if isfield(data, 'time'),       time_val = data.time;
    else,                           time_val = 0.0; 
    end

    % 2. Setup Precision
    if nargin < 3 || isempty(wdsize)
        wdsize = 4; % Default to single for viz speed
        if isa(data.x, 'double'), wdsize = 8; end
    end
    precision = 'single'; if wdsize == 8, precision = 'double'; end

    % 3. Build rdcode
    rdcode = '';
    if isfield(data, 'x'), rdcode = [rdcode 'X']; end
    if isfield(data, 'u'), rdcode = [rdcode 'U']; end
    if isfield(data, 'p'), rdcode = [rdcode 'P']; end
    if isfield(data, 't'), rdcode = [rdcode 'T']; end
    
    npsc = 0;
    while isfield(data, sprintf('s%d', npsc + 1)), npsc = npsc + 1; end
    if npsc > 0, rdcode = [rdcode 'S' num2str(npsc)]; end

    % 4. File I/O
    [path, ~, ~] = fileparts(basepath);
    if ~isempty(path) && ~exist(path, 'dir'), mkdir(path); end

    filename = sprintf('%s0.f%05d', basepath, istep_val + 1);
    [fileID, msg] = fopen(filename, 'wb', 'ieee-le'); 
    assert(isempty(msg), msg);

    % Header (Fixed 132-byte width)
    header = sprintf('#std %1d %2d %2d %2d %10d %10d %20.13E %9d %6d %6d %-10s %15.7E %-22s', ...
        wdsize, nx, ny, nz, nelt, nelt, time_val, istep_val, 0, 1, rdcode, 1.0, 'F');
    header = [header repmat(' ', 1, 132 - length(header))];
    fwrite(fileID, header, 'char');
    fwrite(fileID, 6.54321, 'float32'); % Endian tag

    % 5. Data Blocks
    fwrite(fileID, data.inde, 'int32');

    if isfield(data, 'x')
        xyz = zeros(nxyz, ndim, nelt, precision);
        xyz(:,1,:) = reshape(data.x, [nxyz, 1, nelt]);
        xyz(:,2,:) = reshape(data.y, [nxyz, 1, nelt]);
        if ndim == 3, xyz(:,3,:) = reshape(data.z, [nxyz, 1, nelt]); end
        fwrite(fileID, xyz, precision);
    end

    if isfield(data, 'u')
        uvw = zeros(nxyz, ndim, nelt, precision);
        uvw(:,1,:) = reshape(data.u, [nxyz, 1, nelt]);
        uvw(:,2,:) = reshape(data.v, [nxyz, 1, nelt]);
        if ndim == 3, uvw(:,3,:) = reshape(data.w, [nxyz, 1, nelt]); end
        fwrite(fileID, uvw, precision);
    end

    if isfield(data, 'p'), fwrite(fileID, cast(data.p, precision), precision); end
    if isfield(data, 't'), fwrite(fileID, cast(data.t, precision), precision); end

    for i = 1:npsc
        fwrite(fileID, cast(data.(sprintf('s%d', i)), precision), precision);
    end

    fclose(fileID);
end
