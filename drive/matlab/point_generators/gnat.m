%% --- Optimized Implementation ---
% Fed libROM implementation into Gemini
function [phi] = gnat(Q, m_used, nsr)
    [n, total_m] = size(Q);
    m = min(m_used, total_m);
    ns = iif(nsr > 0, nsr, m);

    phi = zeros(1, ns);
    is_used = false(n, 1);
    base_nsi = floor(ns / m);
    ns_mod_nr = mod(ns, m);
    U = [];
    curr_phi_idx = 0;

    for l = 1:m
        if l == 1
            r = Q(:, 1);
        else
            % Vectorized sampling indexing instead of P matrix multiplication
            sampled_indices = phi(1:curr_phi_idx);
            c = U(sampled_indices, :) \ Q(sampled_indices, l);
            r = Q(:, l) - U * c;
        end
        U = [U, Q(:, l)];
        nsi = base_nsi + (l <= ns_mod_nr);
        for i = 1:nsi
            temp_r = abs(r);
            temp_r(is_used) = -1;
            [~, s_row] = max(temp_r);
            curr_phi_idx = curr_phi_idx + 1;
            phi(curr_phi_idx) = s_row;
            is_used(s_row) = true;
        end
    end
end

function val = iif(cond, t, f), if cond, val = t; else, val = f; end, end
