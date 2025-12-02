function Jacobian = computeJacobian_Boussinesq_optimized(grid, rock, fluid, thermal, physics, ...
                                                connectionsx, connectionsz, ...
                                                y_old, y_new, bc)
% OPTIMIZED VERSION - Uses sparse triplet assembly for better performance
%
% Computes the Jacobian matrix for the Boussinesq convection problem
%
% State vector: [p1; T1; p2; T2; ...; pN; TN]
% Jacobian(i,j) = d(residual_i)/d(y_j)

    N = grid.Nx * grid.Nz;

    % Extract pressure and temperature
    p = y_new(1:2:end);
    T = y_new(2:2:end);

    % Compute density
    rho = fluid.rho_0 * (1 - fluid.beta * (T - physics.T0));

    % Density derivative
    drho_dT = -fluid.rho_0 * fluid.beta;

    % Compute transmissibilities
    [Tx, Tz] = transmissibility(grid, rock, fluid, connectionsx, connectionsz, rho);

    % Compute thermal conductivities
    T_empty = [];
    Ke = thermal.kappa_eff;
    [Kex, Kez] = computeConductivity(grid, rock, fluid, connectionsx, connectionsz, T_empty, Ke);

    cell_volume = grid.dx * grid.dy * grid.dz;

    % Pre-allocate triplet arrays for sparse matrix assembly
    % Estimate: ~20 non-zeros per equation on average (coupling through fluxes)
    max_entries = 20 * 2 * N;
    I = zeros(max_entries, 1);
    J = zeros(max_entries, 1);
    V = zeros(max_entries, 1);
    idx = 0;

    %% ====================================================================
    %  X-DIRECTION CONTRIBUTIONS (VECTORIZED WHERE POSSIBLE)
    %% ====================================================================
    n_connx = size(connectionsx, 1);
    i_cells = connectionsx(:, 1);
    j_cells = connectionsx(:, 2);

    pi_idx = 2*i_cells - 1;
    pj_idx = 2*j_cells - 1;
    Ti_idx = 2*i_cells;
    Tj_idx = 2*j_cells;

    % Pre-compute velocities for all faces
    p_diff = p(j_cells) - p(i_cells);
    face_area = grid.dy * grid.dz;
    u_face_x = -Tx / face_area .* p_diff;

    % Mass equation derivatives (can be added directly to arrays)
    % These are simple and can be vectorized
    new_idx = idx + 1 : idx + 4*n_connx;
    I(new_idx) = [pi_idx; pi_idx; pj_idx; pj_idx];
    J(new_idx) = [pi_idx; pj_idx; pi_idx; pj_idx];
    V(new_idx) = [-Tx; Tx; Tx; -Tx];
    idx = idx + 4*n_connx;

    % Energy equation derivatives (need upwind logic, so loop)
    for f = 1:n_connx
        i = i_cells(f);
        j = j_cells(f);
        pi = pi_idx(f);
        pj = pj_idx(f);
        Ti = Ti_idx(f);
        Tj = Tj_idx(f);

        u_f = u_face_x(f);

        % Determine upwind
        if u_f >= 0
            T_upwind = T(i);
            T_up_idx = Ti;
        else
            T_upwind = T(j);
            T_up_idx = Tj;
        end

        % Advective derivatives
        dAdv_dpi = fluid.rho_0 * fluid.c_f * Tx(f) * T_upwind;
        dAdv_dpj = -fluid.rho_0 * fluid.c_f * Tx(f) * T_upwind;
        dAdv_dTup = fluid.rho_0 * fluid.c_f * u_f * face_area;

        % Conductive derivatives
        dCond_dTi = -Kex(f);
        dCond_dTj =  Kex(f);

        % Add to triplet arrays
        entries = [
            Ti, pi, dAdv_dpi;
            Ti, pj, dAdv_dpj;
            Ti, T_up_idx, dAdv_dTup;
            Ti, Ti, dCond_dTi;
            Ti, Tj, dCond_dTj;
            Tj, pi, -dAdv_dpi;
            Tj, pj, -dAdv_dpj;
            Tj, T_up_idx, -dAdv_dTup;
            Tj, Ti, -dCond_dTi;
            Tj, Tj, -dCond_dTj
        ];

        n_entries = size(entries, 1);
        I(idx+1:idx+n_entries) = entries(:,1);
        J(idx+1:idx+n_entries) = entries(:,2);
        V(idx+1:idx+n_entries) = entries(:,3);
        idx = idx + n_entries;
    end

    %% ====================================================================
    %  Z-DIRECTION CONTRIBUTIONS (with gravity)
    %% ====================================================================
    n_connz = size(connectionsz, 1);
    i_cells = connectionsz(:, 1);
    j_cells = connectionsz(:, 2);

    pi_idx = 2*i_cells - 1;
    pj_idx = 2*j_cells - 1;
    Ti_idx = 2*i_cells;
    Tj_idx = 2*j_cells;

    % Pre-compute for all z-faces
    rho_avg = 0.5 * (rho(i_cells) + rho(j_cells));
    drho_avg_dTi = 0.5 * drho_dT;
    drho_avg_dTj = 0.5 * drho_dT;

    p_diff = p(j_cells) - p(i_cells);
    gravity_term = rho_avg * physics.g * grid.dz;
    face_area = grid.dx * grid.dy;
    u_face_z = -Tz / face_area .* (p_diff - gravity_term);

    % Mass equation derivatives (vectorized)
    % Pressure contributions
    new_idx = idx + 1 : idx + 4*n_connz;
    I(new_idx) = [pi_idx; pi_idx; pj_idx; pj_idx];
    J(new_idx) = [pi_idx; pj_idx; pi_idx; pj_idx];
    V(new_idx) = [-Tz; Tz; Tz; -Tz];
    idx = idx + 4*n_connz;

    % Temperature contributions to mass equation (buoyancy)
    dMass_dTi = -Tz * drho_avg_dTi * physics.g * grid.dz;
    dMass_dTj = -Tz * drho_avg_dTj * physics.g * grid.dz;

    new_idx = idx + 1 : idx + 4*n_connz;
    I(new_idx) = [pi_idx; pi_idx; pj_idx; pj_idx];
    J(new_idx) = [Ti_idx; Tj_idx; Ti_idx; Tj_idx];
    V(new_idx) = [dMass_dTi; dMass_dTj; -dMass_dTi; -dMass_dTj];
    idx = idx + 4*n_connz;

    % Energy equation derivatives (need upwind logic)
    for f = 1:n_connz
        i = i_cells(f);
        j = j_cells(f);
        pi = pi_idx(f);
        pj = pj_idx(f);
        Ti = Ti_idx(f);
        Tj = Tj_idx(f);

        u_f = u_face_z(f);
        rho_avg_f = rho_avg(f);

        % Determine upwind
        if u_f >= 0
            T_upwind = T(i);
            T_up_idx = Ti;
        else
            T_upwind = T(j);
            T_up_idx = Tj;
        end

        % Velocity derivatives
        du_dpi = Tz(f) / face_area;
        du_dpj = -Tz(f) / face_area;
        du_dTi = Tz(f) / face_area * drho_avg_dTi * physics.g * grid.dz;
        du_dTj = Tz(f) / face_area * drho_avg_dTj * physics.g * grid.dz;

        % Advective derivatives
        dAdv_dpi = fluid.rho_0 * fluid.c_f * du_dpi * face_area * T_upwind;
        dAdv_dpj = fluid.rho_0 * fluid.c_f * du_dpj * face_area * T_upwind;
        dAdv_dTup = fluid.rho_0 * fluid.c_f * u_f * face_area;
        dAdv_dTi_vel = fluid.rho_0 * fluid.c_f * du_dTi * face_area * T_upwind;
        dAdv_dTj_vel = fluid.rho_0 * fluid.c_f * du_dTj * face_area * T_upwind;

        % Conductive derivatives
        dCond_dTi = -Kez(f);
        dCond_dTj =  Kez(f);

        % Add to triplet arrays
        entries = [
            Ti, pi, dAdv_dpi;
            Ti, pj, dAdv_dpj;
            Ti, T_up_idx, dAdv_dTup;
            Ti, Ti, dAdv_dTi_vel + dCond_dTi;
            Ti, Tj, dAdv_dTj_vel + dCond_dTj;
            Tj, pi, -dAdv_dpi;
            Tj, pj, -dAdv_dpj;
            Tj, T_up_idx, -dAdv_dTup;
            Tj, Ti, -dAdv_dTi_vel - dCond_dTi;
            Tj, Tj, -dAdv_dTj_vel - dCond_dTj
        ];

        n_entries = size(entries, 1);
        I(idx+1:idx+n_entries) = entries(:,1);
        J(idx+1:idx+n_entries) = entries(:,2);
        V(idx+1:idx+n_entries) = entries(:,3);
        idx = idx + n_entries;
    end

    %% ====================================================================
    %  ACCUMULATION TERM (VECTORIZED)
    %% ====================================================================
    T_cells = (2:2:2*N)';
    dAcc_dT = thermal.c_eff * cell_volume / grid.dt;

    new_idx = idx + 1 : idx + N;
    I(new_idx) = T_cells;
    J(new_idx) = T_cells;
    V(new_idx) = -dAcc_dT;
    idx = idx + N;

    % Trim unused entries and create sparse matrix
    I = I(1:idx);
    J = J(1:idx);
    V = V(1:idx);

    Jacobian = sparse(I, J, V, 2*N, 2*N);

end

