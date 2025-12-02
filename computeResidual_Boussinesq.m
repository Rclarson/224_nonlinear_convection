function residual = computeResidual_Boussinesq(grid, rock, fluid, thermal, physics, ...
                                                connectionsx, connectionsz, ...
                                                y_old, y_new, bc)
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
    T_empty = [];  % Not used in computeConductivity
    Ke = thermal.kappa_eff;
    [Kex, Kez] = computeConductivity(grid, rock, fluid, connectionsx, connectionsz, T_empty, Ke);

    %% ====================================================================
    %  X-DIRECTION FLUXES
    %% ====================================================================
    for f = 1:size(connectionsx, 1)
        i = connectionsx(f, 1);  % Left cell
        j = connectionsx(f, 2);  % Right cell

        % Pressure indices
        pi = 2*i - 1;
        pj = 2*j - 1;

        % Temperature indices
        Ti = 2*i;
        Tj = 2*j;

        % --- Mass conservation (continuity) ---
        % Flux: F = -T * (p_j - p_i)
        mass_flux = Tx(f) * (p(j) - p(i));

        residual(pi) = residual(pi) + mass_flux;
        residual(pj) = residual(pj) - mass_flux;

        % --- Energy equation ---
        % Advective flux (upwind based on velocity direction)
        u_face = -Tx(f) / (grid.dy * grid.dz) * (p(j) - p(i));

        if u_face >= 0  % Flow from i to j
            T_upwind = T(i);
        else  % Flow from j to i
            T_upwind = T(j);
        end

        advective_flux = fluid.rho_0 * fluid.c_f * u_face * grid.dy * grid.dz * T_upwind;

        % Conductive flux
        conductive_flux = Kex(f) * (T(j) - T(i));

        % Total energy flux
        energy_flux = advective_flux + conductive_flux;

        residual(Ti) = residual(Ti) + energy_flux;
        residual(Tj) = residual(Tj) - energy_flux;
    end

    %% ====================================================================
    %  Z-DIRECTION FLUXES (with gravity)
    %% ====================================================================
    for f = 1:size(connectionsz, 1)
        i = connectionsz(f, 1);  % Lower cell
        j = connectionsz(f, 2);  % Upper cell

        % Pressure indices
        pi = 2*i - 1;
        pj = 2*j - 1;

        % Temperature indices
        Ti = 2*i;
        Tj = 2*j;

        % Average density for gravity term
        rho_avg = 0.5 * (rho(i) + rho(j));

        % --- Mass conservation with gravity ---
        % Flux: F = -T * (p_j - p_i - rho_avg * g * dz)
        mass_flux = Tz(f) * (p(j) - p(i) - rho_avg * physics.g * grid.dz);

        residual(pi) = residual(pi) + mass_flux;
        residual(pj) = residual(pj) - mass_flux;

        % --- Energy equation with gravity ---
        % Velocity at face
        u_face = -Tz(f) / (grid.dx * grid.dy) * (p(j) - p(i) - rho_avg * physics.g * grid.dz);

        % Upwind temperature
        if u_face >= 0  % Flow from i to j (upward)
            T_upwind = T(i);
        else  % Flow from j to i (downward)
            T_upwind = T(j);
        end

        advective_flux = fluid.rho_0 * fluid.c_f * u_face * grid.dx * grid.dy * T_upwind;

        % Conductive flux
        conductive_flux = Kez(f) * (T(j) - T(i));

        % Total energy flux
        energy_flux = advective_flux + conductive_flux;

        residual(Ti) = residual(Ti) + energy_flux;
        residual(Tj) = residual(Tj) - energy_flux;
    end

    %% ====================================================================
    %  ACCUMULATION TERM (energy equation only)
    %% ====================================================================
    cell_volume = grid.dx * grid.dy * grid.dz;

    for cell = 1:N
        Ti = 2 * cell;

        % Time derivative: c_eff * V * (T_new - T_old) / dt
        accumulation = thermal.c_eff * cell_volume / grid.dt * (T(cell) - T_old(cell));

        residual(Ti) = residual(Ti) - accumulation;
    end

    % Note: Pressure equations have no accumulation term (incompressible)

end

