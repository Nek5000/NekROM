function test_load_case_metadata()
    % TEST_LOAD_CASE_METADATA Exercises case metadata parsing and fallback.

    driver_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(driver_dir, 'io'));

    shear_dir = fullfile(driver_dir, '..', '..', 'examples', 'shear');

    shear_meta = load_case_metadata(shear_dir, 'shear');
    assert(strcmp(shear_meta.case_file, fullfile(shear_dir, 'shear.mor')), ...
        'Shear case should resolve shear.mor.');
    assert(strcmp(shear_meta.ips, 'H10'), 'Shear case should use H10 metadata from ops/ips.');
    assert(strcmp(shear_meta.pod_mode0, 'avg'), 'Shear case should use avg pod:mode0.');
    assert(shear_meta.subtract_mean, 'Shear case should subtract the mean.');
    assert(strcmp(shear_meta.conv_approach, 'mclsdeim'), ...
        'Shear case should use the mclsdeim runtime mode from shear.mor.');
    assert(shear_meta.nb == 20, 'Shear case should read nb=20.');
    assert(shear_meta.ns == 80, 'Shear case should read ns=80 from ops/ns.');
    assert(shear_meta.deim_nbnl == 20, 'Shear case should read deim:nbnl=20.');
    assert(abs(shear_meta.deim_alpha - 1e-12) < 1e-18, 'Shear case should read deim:alpha=1e-12.');
    assert(shear_meta.deim_dumpnls, 'Shear case should read deim:dumpnls=yes.');
    assert(shear_meta.deim_dumpfine, 'Shear case should read deim:dumpfine=yes.');

    os7000_dir = fullfile(driver_dir, '..', '..', 'examples', 'os7000');
    os7000_meta = load_case_metadata(os7000_dir, 'u3_t020_n13');
    assert(strcmp(os7000_meta.case_file, fullfile(os7000_dir, 'u3_t020_n13.mor')), ...
        'os7000 should resolve u3_t020_n13.mor.');
    assert(strcmp(os7000_meta.ips, 'L2'), 'os7000 should use L2 metadata from the .mor file.');
    assert(strcmp(os7000_meta.pod_mode0, 'state'), 'os7000 should use state pod:mode0.');
    assert(~os7000_meta.subtract_mean, 'State mode0 should disable mean subtraction for os7000.');
    assert(strcmp(os7000_meta.conv_approach, 'mclsdeim'), ...
        'os7000 should use the mclsdeim runtime mode from u3_t020_n13.mor.');
    assert(os7000_meta.nb == 1, 'os7000 should read nb=1.');
    assert(os7000_meta.ns == 100, 'os7000 should read ns=100.');
    assert(os7000_meta.deim_nbnl == 20, 'os7000 should read deim:nbnl=20.');
    assert(abs(os7000_meta.deim_alpha - 1e-12) < 1e-18, 'os7000 should read deim:alpha=1e-12.');
    assert(os7000_meta.deim_dumpnls, 'os7000 should read deim:dumpnls=yes.');
    assert(os7000_meta.deim_dumpfine, 'os7000 should read deim:dumpfine=yes.');

    temp_root = tempname;
    if ~exist(temp_root, 'dir')
        mkdir(temp_root);
    end
    cleanup = onCleanup(@() cleanup_temp_dir(temp_root));

    fallback_dir = fullfile(temp_root, 'fallback_case');
    ops_dir = fullfile(fallback_dir, 'ops');
    if ~exist(fallback_dir, 'dir')
        mkdir(fallback_dir);
    end
    if ~exist(ops_dir, 'dir')
        mkdir(ops_dir);
    end

    mor_text = sprintf([ ...
        '[GENERAL]\n' ...
        'nb = 7\n' ...
        'ns = 13\n' ...
        '\n' ...
        '[POD]\n' ...
        'type = hlm\n' ...
        'mode0 = state\n' ...
        '\n' ...
        '[DEIM]\n' ...
        'mode = clsdeim\n' ...
        'alpha = 2.5e-4\n' ...
        'nbnl = 4\n' ...
        'dumpnls = yes\n' ...
        'dumpfine = no\n']);
    write_text_file(fullfile(fallback_dir, 'bundle.mor'), mor_text);
    write_text_file(fullfile(ops_dir, 'nb'), sprintf('11\n'));
    write_text_file(fullfile(ops_dir, 'ns'), sprintf('19\n'));
    write_text_file(fullfile(ops_dir, 'ips'), sprintf('L2\n'));

    fallback_meta = load_case_metadata(fallback_dir, 'rb_axi');
    assert(strcmp(fallback_meta.case_file, fullfile(fallback_dir, 'bundle.mor')), ...
        'Metadata loader should discover a non-matching .mor filename.');
    assert(fallback_meta.nb == 11, 'ops/nb should override the .mor file.');
    assert(fallback_meta.ns == 19, 'ops/ns should override the .mor file.');
    assert(strcmp(fallback_meta.ips, 'L2'), 'ops/ips should override the .mor file.');
    assert(strcmp(fallback_meta.pod_mode0, 'state'), 'The mode0 field should parse from the .mor file.');
    assert(~fallback_meta.subtract_mean, 'State mode0 should disable mean subtraction.');
    assert(strcmp(fallback_meta.conv_approach, 'clsdeim'), 'deim:mode should parse from the .mor file.');
    assert(abs(fallback_meta.deim_alpha - 2.5e-4) < 1e-18, 'deim:alpha should parse from the .mor file.');
    assert(fallback_meta.deim_nbnl == 4, 'deim:nbnl should parse from the .mor file.');
    assert(fallback_meta.deim_dumpnls, 'deim:dumpnls should parse from the .mor file.');
    assert(~fallback_meta.deim_dumpfine, 'deim:dumpfine should parse from the .mor file.');
end

function write_text_file(path, contents)
    fid = fopen(path, 'w');
    assert(fid >= 0, 'Unable to open %s for writing.', path);
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, '%s', contents);
end

function cleanup_temp_dir(root_dir)
    if exist(root_dir, 'dir')
        rmdir(root_dir, 's');
    end
end
