function save_deim_artifacts(out_dir, rom_data)
    % SAVE_DEIM_ARTIFACTS Persist DEIM-family operators in NekROM text format.

    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    if isfield(rom_data, 'use_full_quadrature') && rom_data.use_full_quadrature
        error('save_deim_artifacts:UnsupportedFullQuadrature', ...
            ['Full-quadrature DEIM data cannot be written to ops/ in the current ' ...
             'Fortran format. Keep this path in memory or save it to a separate diagnostic location.']);
    end

    ndeim_pts = numel(rom_data.inds);
    eval_inds = pick_field(rom_data, 'eval_inds', 'inds_os');
    ndeim_pts_eval = numel(eval_inds);

    if isfield(rom_data, 'inds_os') && ~isempty(rom_data.inds_os)
        inds_os = rom_data.inds_os(:);
    else
        inds_os = eval_inds(:);
    end
    ndeim_pts_os = max(0, numel(inds_os) - ndeim_pts);

    write_int_scalar(fullfile(out_dir, 'deim_npts'), ndeim_pts);
    write_int_scalar(fullfile(out_dir, 'deim_npts_os'), ndeim_pts_os);
    write_int_scalar(fullfile(out_dir, 'deim_npts_eval'), ndeim_pts_eval);

    write_int_vector(fullfile(out_dir, 'deim_inds'), rom_data.inds);
    if ndeim_pts_os > 0
        write_int_vector(fullfile(out_dir, 'deim_inds_os'), inds_os);
    end
    write_int_vector(fullfile(out_dir, 'deim_eval_inds'), eval_inds);

    write_real_vector(fullfile(out_dir, 'deim_eval_weights'), rom_data.eval_weights);

    eval_u_p = pick_field(rom_data, 'eval_u_p', 'u_p');
    eval_v_p = pick_field(rom_data, 'eval_v_p', 'v_p');
    eval_ux_p = pick_field(rom_data, 'eval_ux_p', 'ux_p');
    eval_uy_p = pick_field(rom_data, 'eval_uy_p', 'uy_p');

    write_real_matrix(fullfile(out_dir, 'deim_u_p'), eval_u_p);
    write_real_matrix(fullfile(out_dir, 'deim_v_p'), eval_v_p);
    if isfield(rom_data, 'eval_w_p') || isfield(rom_data, 'w_p')
        eval_w_p = pick_field(rom_data, 'eval_w_p', 'w_p');
        write_real_matrix(fullfile(out_dir, 'deim_w_p'), eval_w_p);
    end
    write_real_matrix(fullfile(out_dir, 'deim_ux_p'), eval_ux_p);
    write_real_matrix(fullfile(out_dir, 'deim_uy_p'), eval_uy_p);
    if isfield(rom_data, 'eval_uz_p') || isfield(rom_data, 'uz_p')
        eval_uz_p = pick_field(rom_data, 'eval_uz_p', 'uz_p');
        write_real_matrix(fullfile(out_dir, 'deim_uz_p'), eval_uz_p);
    end

    write_real_matrix(fullfile(out_dir, 'deim_nl_bas_p_eval'), rom_data.nl_bas_p_eval);
    write_real_matrix(fullfile(out_dir, 'deim_proj_mat'), rom_data.proj_mat);
    write_real_matrix(fullfile(out_dir, 'deim_zmc'), rom_data.zmc);
    write_real_matrix(fullfile(out_dir, 'deim_Ainv'), rom_data.Ainv);
    write_real_matrix(fullfile(out_dir, 'deim_interp_mat'), rom_data.interp_mat);

    if isfield(rom_data, 'mu')
        write_real_vector(fullfile(out_dir, 'deim_mu'), rom_data.mu);
    end
    if isfield(rom_data, 'tau')
        write_real_matrix(fullfile(out_dir, 'deim_tau'), rom_data.tau);
    end
    if isfield(rom_data, 'A_tau_inv')
        write_real_matrix(fullfile(out_dir, 'deim_A_tau_inv'), rom_data.A_tau_inv);
    end
    if isfield(rom_data, 'alpha') && ~isempty(rom_data.alpha)
        write_real_scalar(fullfile(out_dir, 'deim_alpha'), rom_data.alpha);
    end
end

function value = pick_field(struct_data, primary_name, fallback_name)
    if isfield(struct_data, primary_name)
        value = struct_data.(primary_name);
    else
        value = struct_data.(fallback_name);
    end
end

function write_real_vector(path, value)
    fid = fopen(path, 'w');
    assert(fid >= 0, 'Failed to open %s for writing.', path);
    fprintf(fid, '%24.15e\n', value(:));
    fclose(fid);
end

function write_real_matrix(path, value)
    fid = fopen(path, 'w');
    assert(fid >= 0, 'Failed to open %s for writing.', path);
    fprintf(fid, '%24.15e\n', value(:));
    fclose(fid);
end

function write_real_scalar(path, value)
    fid = fopen(path, 'w');
    assert(fid >= 0, 'Failed to open %s for writing.', path);
    fprintf(fid, '%24.15e\n', value);
    fclose(fid);
end

function write_int_scalar(path, value)
    fid = fopen(path, 'w');
    assert(fid >= 0, 'Failed to open %s for writing.', path);
    fprintf(fid, '%d\n', round(value));
    fclose(fid);
end

function write_int_vector(path, value)
    fid = fopen(path, 'w');
    assert(fid >= 0, 'Failed to open %s for writing.', path);
    fprintf(fid, '%d\n', round(value(:)));
    fclose(fid);
end
