function [out_coef] = conv_tensor_sparse(ucoef, pod_u, pod_v, x, y,validator)

    % Create a sparse tensor with arbitary non-zero entries to approximate the
    % convection tensor.
    % Depends on Tensor Toolbox, add it to your MATLABPATH
    % Define a validator function to determine which interactions to include in the
    % sparse tensor.

    persistent tensor

    if isempty(tensor)
        nx1 = size(x,1);
        [zi, w] = zwgll(nx1-1);
        d = deriv_mat(zi);
        [xr,yr,xs,ys,rx,ry,sx,sy,jac,jaci,d] = deriv_geo(x,y,d);
        lgrad=@(u,mode) grad(u,rx,ry,sx,sy,jaci,d,mode);
        Me = reshape(jac.*(w*w'),nL,1);

        nL = prod(size(x));
        nb = size(pod_u,2)

        % Kento's ROM approach. Calculate the gradients of the POD modes
        ux_pods = zeros(size(pod_u));
        uy_pods = zeros(size(pod_u));
        vx_pods = zeros(size(pod_u));
        vy_pods = zeros(size(pod_u));
        
        for i = 1:nb;
            [ux_pod, uy_pod] = lgrad(reshape(pod_u(:,i),size(x)),0);
            [vx_pod, vy_pod] = lgrad(reshape(pod_v(:,i),size(x)),0);
           
            ux_pods(:,i) = reshape(ux_pod, nL,1);
            uy_pods(:,i) = reshape(uy_pod, nL,1);
            vx_pods(:,i) = reshape(vx_pod, nL,1);
            vy_pods(:,i) = reshape(vy_pod, nL,1);
        end;

        % Could also apply Me to pod_u and pod_v inside the loop instead
        % but this seems more efficient.
        pod_weak = [Me.*pod_u(:,2:nb);Me.*pod_v(:,2:nb)];

        % Create the tensor
        tensor = sptensor([nb,nb,nb-1]);
        for i=1:nb;
            for j=1:nb;
                ij = [pod_u(:,i).*ux_pods(:,j) + pod_v(:,i).*uy_pods(:,j); 
                      pod_u(:,i).*vx_pods(:,j) + pod_v(:,i).*vy_pods(:,j)];
                for k=1:nb-1;
                    if validator(i,j,k);
                        tensor(i,j,k) = pod_weak(:,k)'*ij;
                    end;
                end;
            end;
        end
        
    end; 
    % End of if block to persist the tensor. Note, if comparing against the conv_fom approach
    % this will need to be moved back up

    out_coef = double(ttv(tensor, {ucoef, ucoef}, [1,2]));
end
