function run_rom_tests(snaps_path, casename, nb, x_fom, y_fom, pod_u, pod_v, au_full, bu_full, varargin)
% RUN_ROM_TESTS Validates ROM basis and operators against Matlab reference.
%
% Usage:
%   run_rom_tests(path, case, nb, x, y, u, v, au, bu)
%   run_rom_tests(..., 'plot_vort', true)
%   run_rom_tests(..., 'plot_vel_mag', true)

    %% Parse optional plotting arguments
    p = inputParser;
    addParameter(p, 'plot_vort', false, @islogical);
    addParameter(p, 'plot_vel_mag', false, @islogical);
    parse(p, varargin{:});
    
    plot_vort = p.Results.plot_vort;
    plot_vel_mag = p.Results.plot_vel_mag;

    fprintf('--- Starting ROM Validation: %s ---\n', casename);

    %% 1. Generate Reference Data
    % These settings match the original logic for a "ground truth" comparison
    subtract_mean = 1;
    conserve_momentum = 0;
    reorder = 0; 
    method='snapshots';
    inner_product='H10'; % Make sure this matches the NekROM snapshots!

    fprintf('  [1/4] Generating Matlab reference basis...\n');
    snaps_obj = NekSnaps(strcat(snaps_path, casename));
    [pod_ml, u0_full_ml, uk_full_ml] = get_pod_basis(snaps_obj, nb, reorder, subtract_mean, conserve_momentum, method, inner_product);

    % Split reference basis components
    % Assumes pod_ml is a stacked [U; V] matrix
    half_idx = size(pod_ml, 1) / 2;
    pod_u_ml = pod_ml(1:half_idx, 1:nb+1);
    pod_v_ml = pod_ml(half_idx + 1:end, 1:nb+1);

    fprintf('  [2/4] Generating Matlab reference operators (Au, Bu)...\n');
    [au_full_ml, bu_full_ml] = gen_Au(pod_u_ml, pod_v_ml, x_fom, y_fom);

    %% 2. Test: Basis Consistency
    % Compare absolute values to account for arbitrary eigenvector signs
    pod_input = [pod_u; pod_v];
    tol = 1e-5;
    
    err_basis = norm(abs(pod_input) - abs(pod_ml)) / norm(abs(pod_ml));
    assert(err_basis < tol, 'Basis mismatch! Relative error: %e', err_basis);
    fprintf('  [PASS] Basis vectors match (Err: %e)\n', err_basis);

    %% 3. Test: Operator Consistency
    err_bu = norm(abs(bu_full) - abs(bu_full_ml)) / norm(abs(bu_full_ml));
    err_au = norm(abs(au_full) - abs(au_full_ml)) / norm(abs(au_full_ml));

    assert(err_bu < tol, 'Bu operator mismatch! Relative error: %e', err_bu);
    assert(err_au < tol, 'Au operator mismatch! Relative error: %e', err_au);
    fprintf('  [PASS] Operators Au/Bu match (Err Au: %e, Bu: %e)\n', err_au, err_bu);

    %% 4. Reconstruction and Visualization
    if plot_vort || plot_vel_mag
        fprintf('  [4/4] Plotting snapshot reconstructions...\n');
        [u_snaps, v_snaps] = get_snaps(snaps_obj, 1);
        
        fig = figure('Name', ['ROM Reconstruction: ' casename]);
        
        for i = 1:size(uk_full_ml, 2)
            % Current Snapshot Data
            u_proj = u_snaps(:, i);
            v_proj = v_snaps(:, i);

            if plot_vel_mag
                mag = sqrt(u_proj.^2 + v_proj.^2);
                plot_field = reshape(mag, size(x_fom));
                mode_label = 'Velocity Magnitude';
            elseif plot_vort
                % Reshape for spatial derivative calculation
                U_grid = reshape(u_proj, size(x_fom));
                V_grid = reshape(v_proj, size(x_fom));
                plot_field = lcurl(U_grid, V_grid, x_fom, y_fom);
                mode_label = 'Vorticity (\omega)';
            end

            % Visualization
            hold off;
            patch_plot(x_fom, y_fom, plot_field, [], 'PlotType', 'contour');
            title(sprintf('%s - Snapshot %d', mode_label, i));
            colorbar;
            drawnow;
            pause(0.01);
            
            % Check if figure was closed to exit loop gracefully
            if ~ishandle(fig), break; end
        end
    else
        fprintf('  [4/4] Visualization skipped (plots disabled).\n');
    end

    fprintf('--- All tests for %s completed successfully ---\n', casename);
end
