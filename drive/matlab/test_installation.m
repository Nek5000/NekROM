function test_installation()
    % test_installation - Quick sanity check for NekROM MATLAB driver
    %
    % Verifies:
    %   - config.m loads without errors
    %   - All required paths exist
    %   - Key functions are on path
    %   - Octave compatibility (if applicable)
    %
    % Usage:
    %   cd drive/matlab
    %   matlab -batch "test_installation"
    %   octave --eval "test_installation"

    fprintf('\n=== NekROM MATLAB Driver Installation Test ===\n\n');

    driver_dir = fileparts(mfilename('fullpath'));
    cd(driver_dir);

    addpath(fullfile(driver_dir, 'io'));
    addpath(fullfile(driver_dir, 'operators'));
    addpath(fullfile(driver_dir, 'point_generators'));

    toolkit_dir = fullfile(driver_dir, '..', '..', '..', 'NekToolKit', 'matlab');
    if exist(toolkit_dir, 'dir')
        addpath(toolkit_dir, '-begin');
    end

    case_path = '';
    test_count = 0;
    pass_count = 0;

    %% Test 1: Configuration loads
    fprintf('Test 1: Configuration loads... ');
    test_count = test_count + 1;
    try
        config;
        fprintf('[PASS]\n');
        pass_count = pass_count + 1;
    catch ME
        fprintf('[FAIL]\n  Error: %s\n', ME.message);
    end

    %% Test 2: Case metadata loader
    fprintf('Test 2: Case metadata loader... ');
    test_count = test_count + 1;
    try
        test_load_case_metadata();
        fprintf('[PASS]\n');
        pass_count = pass_count + 1;
    catch ME
        fprintf('[FAIL]\n  Error: %s\n', ME.message);
    end

    %% Test 3: Case directory exists
    fprintf('Test 3: Case directory exists (%s)... ', case_path);
    test_count = test_count + 1;
    if exist(case_path, 'dir')
        fprintf('[PASS]\n');
        pass_count = pass_count + 1;
    else
        fprintf('[FAIL]\n  Directory not found: %s\n', case_path);
    end

    %% Test 4: NekToolKit dependency
    fprintf('Test 4: NekToolKit functions available... ');
    test_count = test_count + 1;
    nektoolkit_funcs = {'NekSnaps', 'zwgll', 'deriv_mat', 'deriv_geo', 'grad', 'interp_mat'};
    all_found = true;
    missing = {};
    for i = 1:length(nektoolkit_funcs)
        if exist(nektoolkit_funcs{i}, 'file') ~= 2
            all_found = false;
            missing{end+1} = nektoolkit_funcs{i};
        end
    end
    if all_found
        fprintf('[PASS]\n');
        pass_count = pass_count + 1;
    else
        fprintf('[FAIL]\n  Missing: %s\n', strjoin(missing, ', '));
        fprintf('  Install NekToolKit: https://github.com/kent0/NekToolKit\n');
    end

    %% Test 5: Required modules on path
    fprintf('Test 5: Required modules on path... ');
    test_count = test_count + 1;
    required_functions = {'NekSnaps', 'conv_fom', 'conv_deim', 'qdeim', 'gappy_pod'};
    all_found = true;
    missing = {};
    for i = 1:length(required_functions)
        if exist(required_functions{i}, 'file') ~= 2
            all_found = false;
            missing{end+1} = required_functions{i};
        end
    end
    if all_found
        fprintf('[PASS]\n');
        pass_count = pass_count + 1;
    else
        fprintf('[FAIL]\n  Missing: %s\n', strjoin(missing, ', '));
    end

    %% Test 6: Octave compatibility
    fprintf('Test 6: Octave compatibility check... ');
    test_count = test_count + 1;
    if exist('OCTAVE_VERSION', 'builtin') ~= 0
        % Running in Octave
        try
            pkg load optim;  % Try loading optim package
            fprintf('[PASS] (Octave %s with optim)\n', OCTAVE_VERSION);
            pass_count = pass_count + 1;
        catch
            fprintf('[WARN] (Octave %s without optim)\n', OCTAVE_VERSION);
            fprintf('       Install optim: pkg install -forge optim\n');
            fprintf('       Only needed if ifcopt=true\n');
            pass_count = pass_count + 1;  % Not a hard failure
        end
    else
        % Running in MATLAB
        v = ver('MATLAB');
        fprintf('[PASS] (MATLAB %s)\n', v.Release);
        pass_count = pass_count + 1;
    end

    %% Test 7: Point selection algorithms available
    fprintf('Test 7: Point selection algorithms... ');
    test_count = test_count + 1;
    ps_functions = {'qdeim', 'gappy_pod', 's_opt', 'gnat'};
    ps_found = 0;
    for i = 1:length(ps_functions)
        if exist(ps_functions{i}, 'file') == 2
            ps_found = ps_found + 1;
        end
    end
    if ps_found >= 2  % At least 2 algorithms available
        fprintf('[PASS] (%d/%d available)\n', ps_found, length(ps_functions));
        pass_count = pass_count + 1;
    else
        fprintf('[FAIL] (only %d/%d available)\n', ps_found, length(ps_functions));
    end

    %% Test 8: Archive directory structure
    fprintf('Test 8: Old code archived... ');
    test_count = test_count + 1;
    if exist('archive', 'dir') && ~exist('old', 'dir')
        fprintf('[PASS]\n');
        pass_count = pass_count + 1;
    elseif exist('old', 'dir')
        fprintf('[WARN] "old/" still exists (should be archived)\n');
        pass_count = pass_count + 1;  % Not a hard failure
    else
        fprintf('[PASS] (no deprecated code)\n');
        pass_count = pass_count + 1;
    end

    %% Test 9: README exists
    fprintf('Test 9: Documentation exists... ');
    test_count = test_count + 1;
    if exist('README.md', 'file')
        fprintf('[PASS]\n');
        pass_count = pass_count + 1;
    else
        fprintf('[FAIL] README.md not found\n');
    end

    %% Summary
    fprintf('\n=== Summary ===\n');
    fprintf('Tests passed: %d/%d\n', pass_count, test_count);

    if pass_count == test_count
        fprintf('\n✓ Installation OK - ready to run driver.m\n\n');
        fprintf('Quick start:\n');
        fprintf('  matlab -batch "driver"\n');
        fprintf('  octave --eval "driver"\n\n');
        fprintf('See README.md for full documentation.\n\n');
    elseif pass_count >= test_count - 1
        fprintf('\n⚠ Installation mostly OK - minor issues detected\n\n');
    else
        fprintf('\n✗ Installation has issues - check failures above\n\n');
    end
end
