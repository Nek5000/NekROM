function [bas] = get_pod_basis(snaps_obj, conserve_momentum)
        % Currently only supports the H10 inner product 

        if nargin < 2
            % Off by default
            conserve_momentum = 0
        end

        [u_snaps, v_snaps] = get_snaps(snaps_obj);
        snaps = [u_snaps; v_snaps];
        
        % Calculate mass matrix with Jacobian
        x=snaps_obj.flds{1}.x;
        y=snaps_obj.flds{1}.y;
        nx1 = size(x,1);
        [zi, w] = zwgll(nx1-1);
        d = deriv_mat(zi);
        [xr,yr,xs,ys,rx,ry,sx,sy,jac,jaci,d] = deriv_geo(x,y,d);
        nL = prod(size(x));
        Me = reshape(jac.*(w*w'),nL,1);
        Me_arr = sparse(diag([Me;Me])); 

        if conserve_momentum
        % Momentum conserving basis vectors (for periodic domain at least)
            eu = [ones(size(u_snaps(:,1))); zeros(size(v_snaps(:,1)))];
            ev = [zeros(size(u_snaps(:,1))); ones(size(v_snaps(:,1)))]; 
            E = [eu,ev];
            E = bsxfun(@rdivide, E, sqrt(dot(E, Me_arr*E)));

            % Remove contributions of these basis vectors from snapshots
            snaps = (eye(size(snaps,1)) - E*(E'*Me_arr))*snaps(:,1:2);
        end;

        %% Calculate the POD of the snapshots
        gramian = snaps'*(Me_arr*snaps);
        gramian = 0.5*(gramian + gramian');
        [eigvecs, eigvals] = eig(gramian);
        [eigvals_sorted, sort_inds] = sort(diag(abs(eigvals)),'descend');
        eigvecs = eigvecs(:,sort_inds);
        bas = snaps*eigvecs; % Using all of the modes
        if conserve_momentum; 
            bas = [E,bas];
        end;

        %fileID = fopen('eigvals_sorted.txt', 'w');
        %fprintf(fileID, '%24.15e\n', eigvals_sorted);
        %fclose(fileID);

        % Orthonormalize the basis
        bas = bsxfun(@rdivide, bas, sqrt(dot(bas, Me_arr*bas)));
        
        % Check that the basis is orthonormal in the mass matrix inner product
        %nrm = bas'*(Me_arr*bas)
end
