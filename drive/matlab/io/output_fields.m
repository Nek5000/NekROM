function output_fields(basepath, data, ifvort, ifwrite, ifvis)
    % Pre-fetch size to avoid repeated struct access
    sz = data.size;

    % 1. Vorticity Calculation
    if ifvort
        % Corrected: Passing both u and v components to calculate vorticity
        data.t = lcurl(reshape(data.u, sz), reshape(data.v, sz), data.x, data.y);
    end

    % 2. Visualization
    if ifvis
        hold off;
        if ifvort
            % Use the computed vorticity field
            patch_plot(data.x, data.y, data.t, [], 'PlotType', 'surface');
        else
            % Magnitude calculation: sqrt(u^2 + v^2) is standard for "velocity magnitude"
            % We reshape to the grid size 'sz' for the plot
            vel_mag = reshape(sqrt(data.u.^2 + data.v.^2), sz);
            patch_plot(data.x, data.y, vel_mag, [], 'PlotType', 'surface');
        end
        drawnow; % Forces the plot to update immediately
    end

    % 3. File I/O
    if ifwrite
        fprintf('Writing output %d\n', data.iostep);
        write_field(basepath, data);
    end
end
