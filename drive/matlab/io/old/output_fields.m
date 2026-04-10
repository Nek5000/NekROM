function [] = output_fields(basepath, data, ifvort, ifwrite, ifvis)

    sz = data.size;
    if ifvort;
        data.t = lcurl(reshape(data.u,sz), reshape(data.v,sz), data.x, data.y);%vx-uy;
    end;

    if ifvis; % Set to 1 to enable plotting in MATLAB. This slows the code considerably though.
        hold off;
        % Surface or contour. Note that contours are not supported on the main branch of NekToolkit
        if ifvort
            %patch_plot(data.x,data.y, reshape(data.t,size(x_fom)), [], 'PlotType', 'surface');
            patch_plot(data.x, data.y, data.t, [], 'PlotType', 'surface');
        else
            patch_plot(data.x, data.y, reshape(data.u.^2 + data.v.^2, sz), [], 'PlotType', 'surface');
        end;
    end;
    if ifwrite
        disp(sprintf('Writing output %i', data.iostep));
        write_field(basepath, data);
    end;

end
