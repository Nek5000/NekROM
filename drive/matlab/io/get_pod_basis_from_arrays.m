function [bas, u0, uk] = get_pod_basis(u_snaps, v_snaps, x, y, nb, subtract_mean, conserve_momentum)
        % Currently only supports the H10 inner product 
        % Returns the average as the first column followed by nb basis vectors

        if nargin < 6
            % Off by default
            conserve_momentum = 0;
        end

        nx1 = size(x,1);
        [zi, w] = zwgll(nx1-1);
        d = deriv_mat(zi);
        [xr,yr,xs,ys,rx,ry,sx,sy,jac,jaci,d] = deriv_geo(x,y,d);
        nL = prod(size(x));
        Me = reshape(jac.*(w*w'),nL,1);
        Me_arr = spdiags([Me;Me], 0, 2*nL, 2*nL);
        %Me_arr = sparse(diag([Me;Me]));
        %Me_arr = sparse(diag([Me;Me]*0 + 1));

        % Subtract off the average
        snaps = [u_snaps; v_snaps];
        avg_snaps = mean(snaps, 2);        
        %bnorm_avg = sqrt(avg_snaps'*Me*avg_snaps);
        %avg_snaps = avg_snaps / bnorm_avg;
        %snaps = snaps / bnorm_avg;
        %snaps_orig = snaps;
        if subtract_mean;
            snaps = (snaps - avg_snaps);
        end
        pod_snaps = snaps;
        %snaps = avg_snaps - snaps;

        if conserve_momentum
            % Momentum conserving basis vectors (for periodic domain at least)
            % Only support 2D for the moment.
            eu = [ones(size(u_snaps(:,1))); zeros(size(v_snaps(:,1)))];
            ev = [zeros(size(u_snaps(:,1))); ones(size(v_snaps(:,1)))]; 
            E = [eu,ev];
            E = bsxfun(@rdivide, E, sqrt(dot(E, Me_arr*E)));

            % Remove contributions of these basis vectors from snapshots
            %E
            pod_snaps = (eye(size(snaps,1)) - E*(E'*Me_arr))*snaps;
            %E'*snaps
            %iMe_arr = sparse(diag(1./[Me;Me]));
            %pod_snaps = snaps - iMe_arr*(E*(E'*(Me_arr*snaps)));
            %pod_snaps = snaps - E*(E'*(snaps));
            norm(snaps - pod_snaps)/norm(snaps)
        end;

        %% Calculate the POD of the snapshots
        if 0; % Correlation matrix version
            gramian = pod_snaps'*(Me_arr*pod_snaps);
            gramian = 0.5*(gramian + gramian');
            %[eigvecs, eigvals] = eig(gramian);
            [eigvecs, eigvals] = eigs(gramian, nb, 'largestabs', 'Tolerance', 1e-16);
            % May need to orthogonalize these?
            [eigvals_sorted, sort_inds] = sort(diag(abs(eigvals)),'descend');
            eigvecs = eigvecs(:,sort_inds);
            bas = pod_snaps*eigvecs(:,1:nb); % Using all of the modes
            %bas = orth(bas(:,1:nb), 1e-12);
            bas = bsxfun(@rdivide, bas, sqrt(dot(bas, Me_arr*bas)));
            %exit;
        else % SVD version
            L = sqrt(Me_arr);
            %L = sparse(diag(sqrt([Me;Me])*0 + 1)); 
            [bas,S,V] = svds(L*pod_snaps, nb, 'largest');
            bas = inv(L)*bas;
        end

        if conserve_momentum; 
            bas = [E,bas(:,1:end-2)];
        end;

        %bas'*bas
        %exit

        %fileID = fopen('eigvals_sorted.txt', 'w');
        %fprintf(fileID, '%24.15e\n', eigvals_sorted);
        %fclose(fileID);
        
        % Calculate the coefficients for each time step.
        uk = bas'*Me_arr*snaps; 
        % Check that the basis is orthonormal in the mass matrix inner product
        %nrm = bas'*(Me_arr*bas)

        % Add the average mode to the basis
        if subtract_mean;
            bas = [avg_snaps, bas];
        end;
        uk = [ones(1,size(snaps,2)); uk];
        u0 = uk(:,1); % Hopefully the initial condition is included in the snapshots
end
