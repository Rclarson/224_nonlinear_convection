function residual = computeResidual_Boussinesq_optimized(grid, rock, fluid, thermal, physics, ...
                                                connectionsx, connectionsz, ...
                                                y_old, y_new, bc)
% OPTIMIZED VERSION - Vectorized computations for better performance
%
% Computes the residual for the Boussinesq convection problem
%
% Equations:
%   Mass conservation: div(u) = 0
%   Energy: c_eff*dT/dt + div(rho_0*c_f*u*T) = div(kappa_eff*grad(T))
%
% State vector: [p1; T1; p2; T2; ...; pN; TN]

    N = grid.Nx * grid.Nz;
    residual = zeros(2*N, 1);

    % Extract pressure and temperature
    p = y_new(1:2:end);
    T = y_new(2:2:end);
    T_old = y_old(2:2:end);

    % Compute density (Boussinesq approximation)
    rho = fluid.rho_0 * (1 - fluid.beta * (T - physics.T0));

    % Compute transmissibilities
    [Tx, Tz] = transmissibility(grid, rock, fluid, connectionsx, connectionsz, rho);

    % Compute thermal conductivities
    T_empty = [];
    Ke = thermal.kappa_eff;
    [Kex, Kez] = computeConductivity(grid, rock, fluid, connectionsx, connectionsz, T_empty, Ke);

    %% ====================================================================
    %  X-DIRECTION FLUXES (VECTORIZED)
    %% ====================================================================
    i_cells = connectionsx(:, 1);
    j_cells = connectionsx(:, 2);

    % Pressure indices
    pi_idx = 2*i_cells - 1;
    pj_idx = 2*j_cells - 1;

    % Temperature indices
    Ti_idx = 2*i_cells;
    Tj_idx = 2*j_cells;

    % Mass fluxes (vectorized)
    p_diff = p(j_cells) - p(i_cells);
    mass_flux_x = Tx .* p_diff;

    % Energy fluxes
    face_area = grid.dy * grid.dz;
    u_face_x = -Tx / face_area .* p_diff;

    % Upwind temperature (vectorized)
    T_upwind_x = T(i_cells);
    T_upwind_x(u_face_x < 0) = T(j_cells(u_face_x < 0));

    advective_flux_x = fluid.rho_0 * fluid.c_f * u_face_x * face_area .* T_upwind_x;
    conductive_flux_x = Kex .* (T(j_cells) - T(i_cells));
    energy_flux_x = advective_flux_x + conductive_flux_x;

    % Accumulate fluxes using sparse matrix operations (more efficient)
    for f = 1:length(i_cells)
        residual(pi_idx(f)) = residual(pi_idx(f)) + mass_flux_x(f);
        residual(pj_idx(f)) = residual(pj_idx(f)) - mass_flux_x(f);
        residual(Ti_idx(f)) = residual(Ti_idx(f)) + energy_flux_x(f);
        residual(Tj_idx(f)) = residual(Tj_idx(f)) - energy_flux_x(f);
    end

    %% ====================================================================
    %  Z-DIRECTION FLUXES (VECTORIZED, with gravity)
    %% ====================================================================
    i_cells = connectionsz(:, 1);
    j_cells = connectionsz(:, 2);

    pi_idx = 2*i_cells - 1;
    pj_idx = 2*j_cells - 1;
    Ti_idx = 2*i_cells;
    Tj_idx = 2*j_cells;

    % Average density for gravity term (vectorized)
    rho_avg = 0.5 * (rho(i_cells) + rho(j_cells));

    % Mass fluxes with gravity (vectorized)
    p_diff = p(j_cells) - p(i_cells);
    gravity_term = rho_avg * physics.g * grid.dz;
    mass_flux_z = Tz .* (p_diff - gravity_term);

    % Energy fluxes
    face_area = grid.dx * grid.dy;
    u_face_z = -Tz / face_area .* (p_diff - gravity_term);

    % Upwind temperature (vectorized)
    T_upwind_z = T(i_cells);
    T_upwind_z(u_face_z < 0) = T(j_cells(u_face_z < 0));

    advective_flux_z = fluid.rho_0 * fluid.c_f * u_face_z * face_area .* T_upwind_z;
    conductive_flux_z = Kez .* (T(j_cells) - T(i_cells));
    energy_flux_z = advective_flux_z + conductive_flux_z;

    % Accumulate fluxes
    for f = 1:length(i_cells)
        residual(pi_idx(f)) = residual(pi_idx(f)) + mass_flux_z(f);
        residual(pj_idx(f)) = residual(pj_idx(f)) - mass_flux_z(f);
        residual(Ti_idx(f)) = residual(Ti_idx(f)) + energy_flux_z(f);
        residual(Tj_idx(f)) = residual(Tj_idx(f)) - energy_flux_z(f);
    end

    %% ====================================================================
    %  ACCUMULATION TERM (VECTORIZED)
    %% ====================================================================
    cell_volume = grid.dx * grid.dy * grid.dz;
    T_cells = 2:2:(2*N);

    % Vectorized accumulation
    accumulation = thermal.c_eff * cell_volume / grid.dt * (T - T_old);
    residual(T_cells) = residual(T_cells) - accumulation;

end

