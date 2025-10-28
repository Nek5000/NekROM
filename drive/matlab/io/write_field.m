% See https://nek5000.github.io/NekDoc/problem_setup/case_files.html#restart-output-files-f
% Just handle writing out the coordinates and the velocity field for now.
% This might not be worthwhile since deriv_geo doesn't support 3D currently
% Whatever, might be worth it just to support 2D.

function [] = write_field(basepath, inde, data, sz, time, iostep)

    %TODO: Support 3D

    [path, basename, ~] = fileparts(basepath);
    if prod(size(path)) > 0;
        mkdir(path);
    end;

    % Open file for writing
    filename = sprintf('%s0.f%05d', basepath, iostep + 1);
    [fileID, msg] = fopen(filename, 'W', 'native', 'US-ASCII');
    assert(prod(size(msg)) == 0, msg);

    wdsize = 4; % For visualization, we probably don't need double precision.
    %wdsize = 8;
    if wdsize == 8;
        precision = 'double';
    elseif wdsize == 4;
        precision = 'single';
    else
        print("Invalid wdsize");
        exit;
    end;
    %fldnames = fieldnames(data);
    %first_field = fldnames{1}
    %sz = size(data.(first_field)); 
    % Expand this to support 3D?
    nx = sz(1); 
    ny = sz(2);
    nz = 1;
    nxyz = nx*ny*nz;
    nelt = size(inde,1);
    nelgt = nelt;
    fid = 0;
    nfileoo = 1;
    
    % For now, look at the fields and iostep to see what to write.
    rdcode = '';
    if isfield(data, 'x');
        rdcode='X';
    end;
    if isfield(data, 'u');
        rdcode = append(rdcode, 'U');
    end;
    if isfield(data, 't');
        rdcode = append(rdcode, 'T');
    end;

    ndim = 2;
    
    p0th = 1.0;
    if_press_mesh = 'F';

    % Write the header
    % The specification doesn't mention it, but for visit to read the field files, the
    % numbers need to be right justified and the strings need to be left justified.
    fprintf(fileID, ...
        '#std %1d %2d %2d %2d %10d %10d %20.13E %9d %6d %6d %-10s %15.7E %-22s', ...
        wdsize,nx,ny,nz,nelt,nelgt,time,iostep,fid,nfileoo,rdcode,p0th,if_press_mesh); 
    fwrite(fileID, 6.54321, 'float32');
    assert(ftell(fileID) == 136);

    % Write global element ids
    fwrite(fileID, inde, 'int32');

    % Temporary array to hold data
    tempv = zeros(nxyz,ndim,nelt, precision);

    % Write coordinates
    if contains(rdcode, 'X');
        % Correct
        tempv(:,1,:) = reshape(data.x, [nxyz,nelt]);
        tempv(:,2,:) = reshape(data.y, [nxyz,nelt]);
        fwrite(fileID, tempv, precision);
        % Incorrect
        %fwrite(fileID, x, precision);
        %fwrite(fileID, y, precision);
    end;

    % Write velocity
    if contains(rdcode, 'U');
        tempv(:,1,:) = reshape(data.u, [nxyz,nelt]);
        tempv(:,2,:) = reshape(data.v, [nxyz,nelt]);
        fwrite(fileID, tempv, precision); 
    end;

    % Write pressure
    % TODO

    % Write temperature
    if contains(rdcode, 'T');
        fwrite(fileID, data.t, precision);
        %tempv(:,1,:) = reshape(data.t, [nxyz,nelt]);
        %fwrite(fileID,tempv(:,1,:)
    end;

    % Passive scalars
    % TODO

    % Close file
    fclose(fileID);
end

