function check_dependencies()
% CHECK_DEPENDENCIES Verify all required dependencies for NekROM MATLAB driver
%
% This function checks that:
%   1. Required MATLAB/Octave built-ins are available
%   2. NekToolKit functions are accessible
%   3. Optional features are available (with warnings if missing)
%
% Usage:
%   check_dependencies()  % Run all checks, error if critical deps missing
%
% See also: driver, config

    fprintf('=== NekROM MATLAB Dependency Check ===\n\n');

    % Track overall status
    critical_missing = false;
    optional_missing = false;

    %% Check MATLAB/Octave Environment
    fprintf('Environment:\n');
    if exist('OCTAVE_VERSION', 'builtin')
        fprintf('  ✓ Running on GNU Octave %s\n', OCTAVE_VERSION);
        if octave_package_installed('optim')
            fprintf('  ✓ Octave optim package installed\n');
        else
            fprintf('  ⚠ Octave optim package not found (needed only for ifcopt=true)\n');
            fprintf('    Install with: pkg install -forge optim\n');
            optional_missing = true;
        end
    else
        fprintf('  ✓ Running on MATLAB %s\n', version);
        if exist('fmincon', 'file') == 2
            fprintf('  ✓ fmincon available (Optimization Toolbox)\n');
        else
            fprintf('  ⚠ fmincon not available (needed only for ifcopt=true)\n');
            fprintf('    Install MATLAB Optimization Toolbox if you want constrained optimization.\n');
            optional_missing = true;
        end
    end
    fprintf('\n');

    %% Check NekToolKit (Critical Dependency)
    fprintf('NekToolKit Functions (required):\n');
    required_nektoolkit = {'NekSnaps', 'zwgll', 'deriv_mat', 'deriv_geo', 'grad', 'interp_mat'};

    for i = 1:length(required_nektoolkit)
        fname = required_nektoolkit{i};
        fpath = which(fname);

        if isempty(fpath)
            fprintf('  ✗ MISSING: %s\n', fname);
            critical_missing = true;
        else
            % Check if it's from NekToolKit
            if ~isempty(strfind(fpath, 'NekToolKit'))
                fprintf('  ✓ %s (from NekToolKit)\n', fname);
            else
                fprintf('  ✓ %s (found at: %s)\n', fname, fpath);
            end
        end
    end

    if critical_missing
        fprintf('\n');
        fprintf('ERROR: Missing critical NekToolKit functions.\n');
        fprintf('NekToolKit is required for NekROM MATLAB driver.\n\n');
        fprintf('Installation:\n');
        fprintf('  1. Clone NekToolKit:\n');
        fprintf('     git clone https://github.com/kent0/NekToolKit.git\n');
        fprintf('  2. Add to MATLAB path:\n');
        fprintf('     addpath(''/path/to/NekToolKit/matlab'')\n');
        fprintf('  3. Or place NekToolKit as sibling to NekROM:\n');
        fprintf('     NekROM/\n');
        fprintf('     NekToolKit/\n');
        fprintf('     (driver will auto-detect)\n\n');
        error('check_dependencies:NekToolKitMissing', 'NekToolKit dependency not satisfied.');
    end
    fprintf('\n');

    %% Check Core MATLAB Functions
    fprintf('Core MATLAB/Octave Functions:\n');
    required_core = {'chol', 'svd', 'qr', 'fopen', 'fread', 'fwrite'};

    for i = 1:length(required_core)
        fname = required_core{i};
        if exist(fname, 'builtin') || exist(fname, 'file')
            fprintf('  ✓ %s\n', fname);
        else
            fprintf('  ✗ MISSING: %s\n', fname);
            critical_missing = true;
        end
    end
    fprintf('\n');

    %% Check Optional Features
    fprintf('Optional Features:\n');

    % pagemtimes (MATLAB R2020b+)
    if exist('pagemtimes', 'builtin')
        fprintf('  ✓ pagemtimes (faster tensor assembly)\n');
    else
        fprintf('  ⚠ pagemtimes not available (tensor assembly will be slower)\n');
        fprintf('    Available in MATLAB R2020b+\n');
        optional_missing = true;
    end

    % exportgraphics (MATLAB R2020a+)
    if exist('exportgraphics', 'file')
        fprintf('  ✓ exportgraphics (better PDF plots)\n');
    else
        fprintf('  ⚠ exportgraphics not available (using print fallback)\n');
        fprintf('    Available in MATLAB R2020a+\n');
        optional_missing = true;
    end

    % Parallel Computing Toolbox
    if exist('OCTAVE_VERSION', 'builtin')
        fprintf('  ⚠ parfor not available (Octave does not provide Parallel Computing Toolbox)\n');
        optional_missing = true;
    else
        has_parallel_toolbox = false;
        try
            has_parallel_toolbox = license('test', 'Distrib_Computing_Toolbox') || exist('parpool', 'file') == 2;
        catch
            has_parallel_toolbox = false;
        end

        if has_parallel_toolbox
            fprintf('  ✓ parfor available (Parallel Computing Toolbox)\n');
        else
            fprintf('  ⚠ parfor not available (Parallel Computing Toolbox not installed)\n');
            optional_missing = true;
        end
    end
    fprintf('\n');

    %% Summary
    fprintf('=== Summary ===\n');
    if ~critical_missing && ~optional_missing
        fprintf('✓ All dependencies satisfied!\n');
        fprintf('✓ All optional features available!\n');
    elseif ~critical_missing && optional_missing
        fprintf('✓ All required dependencies satisfied.\n');
        fprintf('⚠ Some optional features unavailable (see above).\n');
        fprintf('  Driver will work but some features may be slower or unavailable.\n');
    else
        fprintf('✗ Critical dependencies missing. Cannot proceed.\n');
    end
    fprintf('\n');
end

function tf = octave_package_installed(pkg_name)
    tf = false;

    try
        packages = pkg('list');
    catch
        return;
    end

    if isempty(packages)
        return;
    end

    if isstruct(packages)
        packages = {packages};
    end

    for i = 1:numel(packages)
        pkg_info = packages{i};
        if isstruct(pkg_info) && isfield(pkg_info, 'name') && strcmp(pkg_info.name, pkg_name)
            tf = true;
            return;
        end
    end
end
