function Jacobian = computeJacobian_Boussinesq(grid, rock, fluid, thermal, physics, ...
                                                connectionsx, connectionsz, ...
                                                y_old, y_new, bc)
% Computes the Jacobian matrix for the Boussinesq convection problem
%
% State vector: [p1; T1; p2; T2; ...; pN; TN]
% Jacobian(i,j) = d(residual_i)/d(y_j)

    N = grid.Nx * grid.Nz;
    Jacobian = sparse(2*N, 2*N);

    % Extract pressure and temperature
    p = y_new(1:2:end);
    T = y_new(2:2:end);

    % Compute density
    rho = fluid.rho_0 * (1 - fluid.beta * (T - physics.T0));

    % Density derivative: drho/dT = -rho_0 * beta
    drho_dT = -fluid.rho_0 * fluid.beta;

    % Compute transmissibilities
    [Tx, Tz] = transmissibility(grid, rock, fluid, connectionsx, connectionsz, rho);

    % Compute thermal conductivities
    T_empty = [];
    Ke = thermal.kappa_eff;
    [Kex, Kez] = computeConductivity(grid, rock, fluid, connectionsx, connectionsz, T_empty, Ke);

    cell_volume = grid.dx * grid.dy * grid.dz;

    %% ====================================================================
    %  X-DIRECTION CONTRIBUTIONS
    %% ====================================================================
    for f = 1:size(connectionsx, 1)
        i = connectionsx(f, 1);
        j = connectionsx(f, 2);

        pi = 2*i - 1;  pj = 2*j - 1;
        Ti = 2*i;      Tj = 2*j;

        % --- Mass equation derivatives ---
        % Mass flux: F = Tx * (p_j - p_i)
        % d(F)/dp_i = -Tx,  d(F)/dp_j = +Tx
        Jacobian(pi, pi) = Jacobian(pi, pi) - Tx(f);
        Jacobian(pi, pj) = Jacobian(pi, pj) + Tx(f);

        Jacobian(pj, pi) = Jacobian(pj, pi) + Tx(f);
        Jacobian(pj, pj) = Jacobian(pj, pj) - Tx(f);

        % --- Energy equation derivatives ---
        % Velocity at face: u = -Tx/(dy*dz) * (p_j - p_i)
        u_face = -Tx(f) / (grid.dy * grid.dz) * (p(j) - p(i));

        % Determine upwind direction
        if u_face >= 0
            T_upwind = T(i);
            T_up_idx = Ti;
        else
            T_upwind = T(j);
            T_up_idx = Tj;
        end

        % Advective flux: F_adv = rho_0 * c_f * u_face * (dy*dz) * T_upwind
        % d(F_adv)/dp_i = rho_0 * c_f * d(u)/dp_i * (dy*dz) * T_upwind
        %               = rho_0 * c_f * (+Tx/(dy*dz)) * (dy*dz) * T_upwind
        %               = rho_0 * c_f * Tx * T_upwind
        % d(F_adv)/dp_j = -rho_0 * c_f * Tx * T_upwind
        dAdv_dpi = fluid.rho_0 * fluid.c_f * Tx(f) * T_upwind;
        dAdv_dpj = -fluid.rho_0 * fluid.c_f * Tx(f) * T_upwind;

        % d(F_adv)/dT_upwind = rho_0 * c_f * u_face * (dy*dz)
        dAdv_dTup = fluid.rho_0 * fluid.c_f * u_face * grid.dy * grid.dz;

        % Conductive flux: F_cond = Kex * (T_j - T_i)
        dCond_dTi = -Kex(f);
        dCond_dTj =  Kex(f);

        % Assemble energy Jacobian (residual += flux)
        Jacobian(Ti, pi) = Jacobian(Ti, pi) + dAdv_dpi;
        Jacobian(Ti, pj) = Jacobian(Ti, pj) + dAdv_dpj;
        Jacobian(Ti, T_up_idx) = Jacobian(Ti, T_up_idx) + dAdv_dTup;
        Jacobian(Ti, Ti) = Jacobian(Ti, Ti) + dCond_dTi;
        Jacobian(Ti, Tj) = Jacobian(Ti, Tj) + dCond_dTj;

        Jacobian(Tj, pi) = Jacobian(Tj, pi) - dAdv_dpi;
        Jacobian(Tj, pj) = Jacobian(Tj, pj) - dAdv_dpj;
        Jacobian(Tj, T_up_idx) = Jacobian(Tj, T_up_idx) - dAdv_dTup;
        Jacobian(Tj, Ti) = Jacobian(Tj, Ti) - dCond_dTi;
        Jacobian(Tj, Tj) = Jacobian(Tj, Tj) - dCond_dTj;
    end

    %% ====================================================================
    %  Z-DIRECTION CONTRIBUTIONS (with gravity)
    %% ====================================================================
    for f = 1:size(connectionsz, 1)
        i = connectionsz(f, 1);  % Lower cell
        j = connectionsz(f, 2);  % Upper cell

        pi = 2*i - 1;  pj = 2*j - 1;
        Ti = 2*i;      Tj = 2*j;

        rho_avg = 0.5 * (rho(i) + rho(j));
        drho_avg_dTi = 0.5 * drho_dT;
        drho_avg_dTj = 0.5 * drho_dT;

        % --- Mass equation derivatives with gravity ---
        % Mass flux: F_mass = Tz * (p_j - p_i - rho_avg * g * dz)

        % Derivatives w.r.t. pressure
        % d(F_mass)/dp_i = -Tz
        % d(F_mass)/dp_j = +Tz
        Jacobian(pi, pi) = Jacobian(pi, pi) - Tz(f);
        Jacobian(pi, pj) = Jacobian(pi, pj) + Tz(f);

        Jacobian(pj, pi) = Jacobian(pj, pi) + Tz(f);
        Jacobian(pj, pj) = Jacobian(pj, pj) - Tz(f);

        % Derivatives w.r.t. temperature (through buoyancy: -rho_avg*g*dz)
        % d(F_mass)/dT_i = Tz * (-drho_avg/dT_i * g * dz) = -Tz * drho_avg_dTi * g * dz
        % d(F_mass)/dT_j = Tz * (-drho_avg/dT_j * g * dz) = -Tz * drho_avg_dTj * g * dz
        dMass_dTi = -Tz(f) * drho_avg_dTi * physics.g * grid.dz;
        dMass_dTj = -Tz(f) * drho_avg_dTj * physics.g * grid.dz;

        Jacobian(pi, Ti) = Jacobian(pi, Ti) + dMass_dTi;
        Jacobian(pi, Tj) = Jacobian(pi, Tj) + dMass_dTj;

        Jacobian(pj, Ti) = Jacobian(pj, Ti) - dMass_dTi;
        Jacobian(pj, Tj) = Jacobian(pj, Tj) - dMass_dTj;

        % --- Energy equation derivatives ---
        % Velocity at face: u_z = -Tz/(dx*dy) * (p_j - p_i - rho_avg*g*dz)
        u_face = -Tz(f) / (grid.dx * grid.dy) * (p(j) - p(i) - rho_avg * physics.g * grid.dz);

        % Determine upwind direction
        if u_face >= 0
            T_upwind = T(i);
            T_up_idx = Ti;
        else
            T_upwind = T(j);
            T_up_idx = Tj;
        end

        % Derivatives of velocity
        % d(u_z)/dp_i = -Tz/(dx*dy) * (-1) = +Tz/(dx*dy)
        % d(u_z)/dp_j = -Tz/(dx*dy) * (+1) = -Tz/(dx*dy)
        du_dpi = Tz(f) / (grid.dx * grid.dy);
        du_dpj = -Tz(f) / (grid.dx * grid.dy);

        % d(u_z)/dT_i = -Tz/(dx*dy) * (-drho_avg_dTi * g * dz)
        du_dTi = Tz(f) / (grid.dx * grid.dy) * drho_avg_dTi * physics.g * grid.dz;
        du_dTj = Tz(f) / (grid.dx * grid.dy) * drho_avg_dTj * physics.g * grid.dz;

        % Advective flux: F_adv = rho_0 * c_f * u_z * (dx*dy) * T_upwind
        % d(F_adv)/dp_i = rho_0 * c_f * d(u_z)/dp_i * (dx*dy) * T_upwind
        dAdv_dpi = fluid.rho_0 * fluid.c_f * du_dpi * grid.dx * grid.dy * T_upwind;
        dAdv_dpj = fluid.rho_0 * fluid.c_f * du_dpj * grid.dx * grid.dy * T_upwind;

        % d(F_adv)/dT_upwind = rho_0 * c_f * u_z * (dx*dy)
        dAdv_dTup = fluid.rho_0 * fluid.c_f * u_face * grid.dx * grid.dy;

        % d(F_adv)/dT_i (from velocity dependence on temperature through buoyancy)
        dAdv_dTi_vel = fluid.rho_0 * fluid.c_f * du_dTi * grid.dx * grid.dy * T_upwind;
        dAdv_dTj_vel = fluid.rho_0 * fluid.c_f * du_dTj * grid.dx * grid.dy * T_upwind;

        % Conductive flux: F_cond = Kez * (T_j - T_i)
        dCond_dTi = -Kez(f);
        dCond_dTj =  Kez(f);

        % Assemble energy Jacobian
        Jacobian(Ti, pi) = Jacobian(Ti, pi) + dAdv_dpi;
        Jacobian(Ti, pj) = Jacobian(Ti, pj) + dAdv_dpj;
        Jacobian(Ti, T_up_idx) = Jacobian(Ti, T_up_idx) + dAdv_dTup;
        Jacobian(Ti, Ti) = Jacobian(Ti, Ti) + dAdv_dTi_vel + dCond_dTi;
        Jacobian(Ti, Tj) = Jacobian(Ti, Tj) + dAdv_dTj_vel + dCond_dTj;

        Jacobian(Tj, pi) = Jacobian(Tj, pi) - dAdv_dpi;
        Jacobian(Tj, pj) = Jacobian(Tj, pj) - dAdv_dpj;
        Jacobian(Tj, T_up_idx) = Jacobian(Tj, T_up_idx) - dAdv_dTup;
        Jacobian(Tj, Ti) = Jacobian(Tj, Ti) - dAdv_dTi_vel - dCond_dTi;
        Jacobian(Tj, Tj) = Jacobian(Tj, Tj) - dAdv_dTj_vel - dCond_dTj;
    end

    %% ====================================================================
    %  ACCUMULATION TERM (energy equation only)
    %% ====================================================================
    for cell = 1:N
        Ti = 2 * cell;

        % d(accumulation)/dT = c_eff * V / dt
        dAcc_dT = thermal.c_eff * cell_volume / grid.dt;

        Jacobian(Ti, Ti) = Jacobian(Ti, Ti) - dAcc_dT;
    end

end

