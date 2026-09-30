function [sort_order] = get_sort_order(x, y)

    indices = [];
    phys_co = []; % Physical coordinates of elements
    for ie = 1:size(x,3)
        % Not necessarily super robust, but should work for now.
        phys_co = [phys_co;[ie, mean(x(:,:,ie), "all"), mean(y(:,:,ie), "all")]];
    end; 
    [~, I] = sort(phys_co(:,2));
    phys_co = phys_co(I,:);
    [~, I] = sort(phys_co(:,3));
    phys_co = phys_co(I,:);
    sort_order = phys_co(:,1);
end
