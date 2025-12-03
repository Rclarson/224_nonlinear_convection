% Boussinesq Convection in 2D Porous Medium
% Coupled pressure-temperature system with buoyancy-driven thermal
% convection
%
% Governing equations:
%   u = -K/mu * (grad(p) - rho_0*(1 - beta*(T-T0))*g)
%   div(u) = 0
%   c_eff*dT/dt + div(rho_0*c_f*u*T) - div(kappa_eff*grad(T)) + Q = 0
%
% Boundary conditions:
%   T = T_h at bottom, T = T_c at top
%   dT/dn = 0 on sides (thermal insulated)
%   u*n = 0 on all walls (no flow)

clc; clear; close all;

%% =======================================================================
%   SIMULATION PARAMETERS
%% =======================================================================
% This structure contains all the parameters for the simulation
% Modify these values to change the problem setup

problem_data = struct();

% --- DOMAIN GEOMETRY ---
problem_data.domain_width  = 200.0;        % Domain width [m]
problem_data.domain_height = 100.0;        % Domain height [m]

% --- GRID RESOLUTION ---
problem_data.n_cells_x = 200;            % Number of cells in x-direction
problem_data.n_cells_z = 100;             % Number of cells in z-direction

% --- ROCK PROPERTIES ---
problem_data.permeability_x = 1e-12;     % Permeability in x [m^2]
problem_data.permeability_z = 1e-12;     % Permeability in z [m^2]
problem_data.porosity       = 0.3;       % Porosity [-]

% --- FLUID PROPERTIES ---
problem_data.fluid_viscosity     = 1e-3;    % Dynamic viscosity [Pa·s]
problem_data.fluid_density_ref   = 1000;    % Reference density [kg/m^3]
problem_data.thermal_expansion   = 2e-4;    % Thermal expansion coefficient [1/K]
problem_data.fluid_heat_capacity = 4180;    % Specific heat capacity [J/(kg·K)]

% --- THERMAL PROPERTIES ---
problem_data.thermal_conductivity = 2.0;    % Effective thermal conductivity [W/(m·K)]

% --- PHYSICS ---
problem_data.gravity            = 9.81;     % Gravitational acceleration [m/s^2]
problem_data.reference_temp     = 300;      % Reference temperature [K]

% --- BOUNDARY CONDITIONS ---
%Conductive solution 
%problem_data.temp_bottom = 350;             % Bottom temperature (hot) [K] (low: Ra -> 2.05e+01)

%Onset of the thermal convection
%problem_data.temp_bottom = 400;             % Bottom temperature (hot) [K] (onset: Ra -> 4.10e+01 )

%Moderate thermal convection
problem_data.temp_bottom = 425;             % Bottom temperature (hot) [K] (mod: Ra -> 5.13e+01)

%High thermal convection
%problem_data.temp_bottom = 450;             % Bottom temperature (hot) [K] (high: Ra -> 6.15e+01)

problem_data.temp_top    = 300;             % Top temperature (cold) [K]
problem_data.pressure_ref = 1e5;            % Reference pressure [Pa]

% --- INITIAL CONDITIONS ---
problem_data.perturbation_amplitude = 0.001;    % Temperature perturbation amplitude [K]
problem_data.perturbation_modes     = 2;        % Number of sinusoidal modes in x

% --- TIME INTEGRATION ---
day = 86400;
year = 365 * day;
problem_data.final_time = 1000 * year;             % Final simulation time [s]
problem_data.time_step  = 10 * year;              % Time step size [s]

% --- SOLVER PARAMETERS ---
problem_data.newton_max_iter = 15;          % Maximum Newton iterations
problem_data.newton_tol_du   = 1e-8;        % Tolerance on update norm
problem_data.newton_tol_res  = 1e-6;        % Tolerance on residual norm

% --- DIAGNOSTICS ---
problem_data.print_interval = 10;           % Print diagnostics every N steps

% --- DIMENSIONLESS NUMBERS ---
% Calculate Rayleigh number for porous media convection
% Ra = (rho_0 * beta * g * K * Delta_T * H) / (mu * alpha)
% where alpha = kappa_eff / (rho_0 * c_f) is thermal diffusivity
Delta_T = problem_data.temp_bottom - problem_data.temp_top;
alpha = problem_data.thermal_conductivity / ...
        (problem_data.fluid_density_ref * problem_data.fluid_heat_capacity);

problem_data.rayleigh_number = ...
    (problem_data.fluid_density_ref * problem_data.thermal_expansion * ...
     problem_data.gravity * problem_data.permeability_z * Delta_T * ...
     problem_data.domain_height) / (problem_data.fluid_viscosity * alpha);

problem_data.thermal_diffusivity = alpha;   % Thermal diffusivity [m^2/s]

% Create output folder
output_folder = 'output_folder';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
fprintf('Output folder created: %s\n\n', output_folder);

grid.Lx = problem_data.domain_width;
grid.Lz = problem_data.domain_height;
grid.Ly = 1.0;  % Not used in 2D, kept for compatibility

grid.Nx = problem_data.n_cells_x;
grid.Nz = problem_data.n_cells_z;
grid.Ny = 1;    % Always 1 for 2D

% Derived grid quantities
grid.dx = grid.Lx / grid.Nx;
grid.dz = grid.Lz / grid.Nz;
grid.dy = grid.Ly / grid.Ny;

N = grid.Nx * grid.Nz;  % Total number of cells

rock.kx  = problem_data.permeability_x * ones(N, 1);
rock.kz  = problem_data.permeability_z * ones(N, 1);
rock.phi = problem_data.porosity * ones(N, 1);

fluid.mu    = problem_data.fluid_viscosity;
fluid.rho_0 = problem_data.fluid_density_ref;
fluid.beta  = problem_data.thermal_expansion;
fluid.c_f   = problem_data.fluid_heat_capacity;

thermal.kappa_eff = problem_data.thermal_conductivity;
thermal.c_eff     = fluid.rho_0 * fluid.c_f * rock.phi(1);

physics.g  = problem_data.gravity;
physics.T0 = problem_data.reference_temp;

bc.T_hot  = problem_data.temp_bottom;
bc.T_cold = problem_data.temp_top;
bc.p_ref  = problem_data.pressure_ref;

% Compute bottom and top cell indices
bc.bottom_cells = 1:grid.Nx;
bc.top_cells    = (grid.Nz-1)*grid.Nx + (1:grid.Nx);

% Reference pressure cell (middle of domain)
bc.p_ref_cell = round(N/2);

% Compute cell coordinates
z_coords = zeros(N, 1);
x_coords = zeros(N, 1);
for j = 1:grid.Nz
    for i = 1:grid.Nx
        cell_idx = (j-1)*grid.Nx + i;
        x_coords(cell_idx) = (i - 0.5) * grid.dx;
        z_coords(cell_idx) = (j - 0.5) * grid.dz;
    end
end

% Linear temperature profile from top to bottom
T_init = bc.T_cold + (bc.T_hot - bc.T_cold) * (1 - z_coords/grid.Lz);

% sinusoidal perturbation to trigger convection
T_init = T_init + problem_data.perturbation_amplitude * ...
         sin(2 * pi * problem_data.perturbation_modes * x_coords / grid.Lx);

% Initial pressure (hydrostatic reference)
p_init = bc.p_ref * ones(N, 1);

% Combine into state vector: [p1; T1; p2; T2; ...; pN; TN]
y_init = zeros(2*N, 1);
y_init(1:2:end) = p_init;
y_init(2:2:end) = T_init;


time.t_final = problem_data.final_time;
time.dt      = problem_data.time_step;
time.t       = 0;

grid.dt    = time.dt;
grid.t_end = time.t_final;


solver.max_iter = problem_data.newton_max_iter;
solver.tol_du   = problem_data.newton_tol_du;
solver.tol_res  = problem_data.newton_tol_res;

connectionsx = connectionsListx(grid);
connectionsz = connectionsListz(grid);

% Calculate simulation parameters
n_time_steps = ceil(time.t_final / time.dt);
n_dof = 2 * N;  % 2 unknowns per cell (p, T)

fprintf('=======================================================\n');
fprintf('  2D Boussinesq Convection Simulation\n');
fprintf('=======================================================\n');
fprintf('Domain: %.2f x %.2f m\n', grid.Lx, grid.Lz);
fprintf('Grid:   %d x %d cells\n', grid.Nx, grid.Nz);
fprintf('Total cells: %d\n', N);
fprintf('Degrees of freedom: %d\n', n_dof);
fprintf('Delta T: %.2f K\n', bc.T_hot - bc.T_cold);
fprintf('Rayleigh number: %.2e\n', problem_data.rayleigh_number);
fprintf('Time step size: %.2f years\n', time.dt / year);
fprintf('Final time: %.2f years\n', time.t_final / year);
fprintf('Expected time steps: %d\n', n_time_steps);
fprintf('Output folder: %s\n', output_folder);
fprintf('=======================================================\n\n');

y_current = y_init;
y_old     = y_init;

step = 1;

% Optimize memory by storing only snapshots instead of full history
n_snapshots = min(10, n_time_steps + 1);
snapshot_interval = max(1, floor(n_time_steps / (n_snapshots - 1)));
t_steps = zeros(1, n_snapshots);
y_history = zeros(2*N, n_snapshots);
snapshot_idx = 1;

% Store initial condition
t_steps(snapshot_idx) = time.t;
y_history(:, snapshot_idx) = y_current;
snapshot_idx = snapshot_idx + 1;

while time.t < time.t_final

    fprintf('Time step %d: t = %.2f years\n', step, time.t / year);

    % Newton iteration
    y_new = y_current;
    converged = false;

    for iter = 1:solver.max_iter

        % Compute residual (using optimized version)
        residual = computeResidual_Boussinesq_optimized(grid, rock, fluid, thermal, physics, ...
                                              connectionsx, connectionsz, ...
                                              y_old, y_new, bc);

        % Compute Jacobian (using optimized version)
        Jacobian = computeJacobian_Boussinesq_optimized(grid, rock, fluid, thermal, physics, ...
                                              connectionsx, connectionsz, ...
                                              y_old, y_new, bc);

        % Apply boundary conditions
        [residual, Jacobian] = applyBC(residual, Jacobian, y_new, grid, bc);

        % Solve for Newton update
        dy = -Jacobian \ residual;

        % Update solution
        y_new = y_new + dy;

        % Check convergence
        norm_du = norm(dy);
        norm_res = norm(residual);

        if norm_du < solver.tol_du && norm_res < solver.tol_res
            converged = true;
            fprintf('  Iter %2d: ||dx|| = %.3e, ||res|| = %.3e\n', ...
                    iter, norm_du, norm_res);
            fprintf('  --> Converged in %d iterations\n', iter);
            break;
        end
    end

    if ~converged
        warning('Newton solver did not converge at t = %.2f', time.t);
    end

    % Accept time step
    time.t = time.t + time.dt;
    y_old = y_new;
    y_current = y_new;

    step = step + 1;

    % Store snapshots at regular intervals
    if mod(step - 1, snapshot_interval) == 0 && snapshot_idx <= n_snapshots
        t_steps(snapshot_idx) = time.t;
        y_history(:, snapshot_idx) = y_current;
        snapshot_idx = snapshot_idx + 1;
    end

    % Compute diagnostics every N steps
    if mod(step, problem_data.print_interval) == 0
        T_map = y_current(2:2:end);
        T_mean = mean(T_map);
        T_min = min(T_map);
        T_max = max(T_map);
        fprintf('  T range: [%.2f, %.2f] K, mean = %.2f K (t = %.2f years)\n', ...
                T_min, T_max, T_mean, time.t / year);
    end
end

% Store final state if not already stored
if snapshot_idx <= n_snapshots
    t_steps(snapshot_idx) = time.t;
    y_history(:, snapshot_idx) = y_current;
    snapshot_idx = snapshot_idx + 1;
end

% Trim unused snapshot slots
if snapshot_idx <= n_snapshots
    t_steps = t_steps(1:snapshot_idx-1);
    y_history = y_history(:, 1:snapshot_idx-1);
end

fprintf('\n=======================================================\n');
fprintf('  Simulation completed\n');
fprintf('=======================================================\n');
fprintf('Total time steps completed: %d\n', step);
fprintf('Final simulation time: %.2f years\n', time.t / year);
fprintf('=======================================================\n\n');

%% =======================================================================
%  POST-PROCESSING AND VISUALIZATION
%% =======================================================================

% Extract final pressure and temperature fields
p_final = y_history(1:2:end, end);
T_final = y_history(2:2:end, end);
T_initial = y_history(2:2:end, 1);

% Compute temperature difference
Delta_T = T_final - T_initial;

% Reshape to 2D grids for plotting
P_map = reshape(p_final, grid.Nx, grid.Nz)';
T_map = reshape(T_final, grid.Nx, grid.Nz)';
Delta_T_map = reshape(Delta_T, grid.Nx, grid.Nz)';

% Create coordinate arrays
x = linspace(grid.dx/2, grid.Lx - grid.dx/2, grid.Nx);
z = linspace(grid.dz/2, grid.Lz - grid.dz/2, grid.Nz);
[X, Z] = meshgrid(x, z);

% Compute velocity field
[u_x, u_z] = computeVelocity(y_history(:,end), grid, rock, fluid, physics, ...
                              connectionsx, connectionsz);

U_map = reshape(u_x, grid.Nx, grid.Nz)';
W_map = reshape(u_z, grid.Nx, grid.Nz)';

%% Plot pressure evolution (separate plots)
fprintf('\nGenerating pressure evolution plots...\n');
n_plots = size(y_history, 2);
time_indices = 1:n_plots;

for i = 1:n_plots
    fig_press = figure('Position', [100, 100, 800, 600]);
    p_snap = reshape(y_history(1:2:end, time_indices(i)), grid.Nx, grid.Nz)';
    contourf(X, Z, flipud(p_snap), 20, 'LineColor', 'none');
    colorbar;
    colormap(jet);
    xlabel('x [m]', 'FontSize', 12);
    ylabel('z [m]', 'FontSize', 12);
    title(sprintf('Pressure at t = %.2f years', t_steps(time_indices(i)) / year), 'FontSize', 14);
    axis equal tight;
    set(gca, 'FontSize', 11);

    % Save with time suffix
    time_suffix = sprintf('t_%05.0f', t_steps(time_indices(i)));
    saveas(fig_press, fullfile(output_folder, ['pressure_evolution_' time_suffix '.png']));
    fprintf('  Saved: pressure_evolution_%s.png\n', time_suffix);

    close(fig_press);
end

%% Plot temperature evolution (separate plots)
fprintf('\nGenerating temperature evolution plots...\n');

for i = 1:n_plots
    fig_temp = figure('Position', [100, 600, 800, 600]);
    T_snap = reshape(y_history(2:2:end, time_indices(i)), grid.Nx, grid.Nz)';
    contourf(X, Z, flipud(T_snap), 20, 'LineColor', 'none');
    colorbar;
    colormap(jet);
    caxis([bc.T_cold, bc.T_hot]);
    xlabel('x [m]', 'FontSize', 12);
    ylabel('z [m]', 'FontSize', 12);
    title(sprintf('Temperature at t = %.2f years', t_steps(time_indices(i)) / year), 'FontSize', 14);
    axis equal tight;
    set(gca, 'FontSize', 11);

    % Save with time suffix
    time_suffix = sprintf('t_%05.0f', t_steps(time_indices(i)));
    saveas(fig_temp, fullfile(output_folder, ['temperature_evolution_' time_suffix '.png']));
    fprintf('  Saved: temperature_evolution_%s.png\n', time_suffix);

    close(fig_temp);
end

%% Plot temperature difference evolution (separate plots)
fprintf('\nGenerating temperature difference evolution plots...\n');

for i = 1:n_plots
    fig_delta = figure('Position', [100, 1100, 800, 600]);
    T_snap = reshape(y_history(2:2:end, time_indices(i)), grid.Nx, grid.Nz)';
    T_init_map = reshape(y_history(2:2:end, 1), grid.Nx, grid.Nz)';
    Delta_T_snap = T_snap - T_init_map;
    contourf(X, Z, flipud(Delta_T_snap), 20, 'LineColor', 'none');
    colorbar;
    colormap(jet);
    xlabel('x [m]', 'FontSize', 12);
    ylabel('z [m]', 'FontSize', 12);
    title(sprintf('\\Delta T = T - T_{init} at t = %.2f years', t_steps(time_indices(i)) / year), 'FontSize', 14);
    axis equal tight;
    set(gca, 'FontSize', 11);

    % Save with time suffix
    time_suffix = sprintf('t_%05.0f', t_steps(time_indices(i)));
    saveas(fig_delta, fullfile(output_folder, ['deltaT_evolution_' time_suffix '.png']));
    fprintf('  Saved: deltaT_evolution_%s.png\n', time_suffix);

    close(fig_delta);
end

%% Plot velocity fields at different times (separate plots)
fprintf('\nGenerating velocity field plots...\n');
skip = 3;  % Downsample factor for velocity arrows

for i = 1:n_plots
    % Extract temperature and compute velocity for this snapshot
    y_snap = y_history(:, time_indices(i));
    T_snap = reshape(y_snap(2:2:end), grid.Nx, grid.Nz)';

    % Compute velocity field at this time
    [u_x_snap, u_z_snap] = computeVelocity(y_snap, grid, rock, fluid, physics, ...
                                           connectionsx, connectionsz);
    U_snap = reshape(u_x_snap, grid.Nx, grid.Nz)';
    W_snap = reshape(u_z_snap, grid.Nx, grid.Nz)';

    % Create individual figure for this time step
    fig_vel = figure('Position', [100, 100, 800, 600]);
    contourf(X, Z, flipud(T_snap), 20, 'LineColor', 'none');
    hold on;
    quiver(X(1:skip:end, 1:skip:end), Z(1:skip:end, 1:skip:end), ...
           U_snap(1:skip:end, 1:skip:end), flipud(W_snap(1:skip:end, 1:skip:end)), ...
           2, 'k', 'LineWidth', 1.2);
    colorbar;
    colormap(jet);
    caxis([bc.T_cold, bc.T_hot]);
    xlabel('x [m]', 'FontSize', 12);
    ylabel('z [m]', 'FontSize', 12);
    title(sprintf('Temperature + Velocity Field at t = %.2f years', t_steps(time_indices(i)) / year), 'FontSize', 14);
    axis equal tight;
    set(gca, 'FontSize', 11);

    % Save with time suffix
    time_suffix = sprintf('t_%05.0f', t_steps(time_indices(i)));
    saveas(fig_vel, fullfile(output_folder, ['velocity_field_' time_suffix '.png']));
    fprintf('  Saved: velocity_field_%s.png\n', time_suffix);

    close(fig_vel);
end

%% Save problem_data structure
save(fullfile(output_folder, 'problem_data.mat'), 'problem_data');
fprintf('Saved: problem_data.mat\n');

%% Save simulation summary report
report_file = fullfile(output_folder, 'simulation_report.txt');
fid = fopen(report_file, 'w');

fprintf(fid, '=======================================================\n');
fprintf(fid, '  2D Boussinesq Convection Simulation Report\n');
fprintf(fid, '=======================================================\n');
fprintf(fid, 'Date: %s\n\n', datestr(now));

fprintf(fid, 'DOMAIN AND GRID:\n');
fprintf(fid, '  Domain size: %.2f x %.2f m\n', problem_data.domain_width, problem_data.domain_height);
fprintf(fid, '  Grid resolution: %d x %d cells\n', problem_data.n_cells_x, problem_data.n_cells_z);
fprintf(fid, '  Total number of cells: %d\n', N);
fprintf(fid, '  Cell dimensions: dx = %.4f m, dz = %.4f m\n', grid.dx, grid.dz);
fprintf(fid, '  Degrees of freedom: %d\n\n', n_dof);

fprintf(fid, 'PHYSICAL PARAMETERS:\n');
fprintf(fid, '  Permeability (x): %.2e m^2\n', problem_data.permeability_x);
fprintf(fid, '  Permeability (z): %.2e m^2\n', problem_data.permeability_z);
fprintf(fid, '  Porosity: %.2f\n', problem_data.porosity);
fprintf(fid, '  Fluid viscosity: %.2e Pa·s\n', problem_data.fluid_viscosity);
fprintf(fid, '  Reference density: %.1f kg/m^3\n', problem_data.fluid_density_ref);
fprintf(fid, '  Thermal expansion: %.2e 1/K\n', problem_data.thermal_expansion);
fprintf(fid, '  Fluid heat capacity: %.1f J/(kg·K)\n', problem_data.fluid_heat_capacity);
fprintf(fid, '  Effective conductivity: %.2f W/(m·K)\n', problem_data.thermal_conductivity);
fprintf(fid, '  Effective heat capacity: %.2e J/(m^3·K)\n', thermal.c_eff);
fprintf(fid, '  Gravity: %.2f m/s^2\n', problem_data.gravity);
fprintf(fid, '  Reference temperature: %.2f K\n\n', problem_data.reference_temp);

fprintf(fid, 'BOUNDARY CONDITIONS:\n');
fprintf(fid, '  Bottom temperature: %.2f K\n', problem_data.temp_bottom);
fprintf(fid, '  Top temperature: %.2f K\n', problem_data.temp_top);
fprintf(fid, '  Temperature difference: %.2f K\n', problem_data.temp_bottom - problem_data.temp_top);
fprintf(fid, '  Reference pressure: %.2e Pa\n\n', problem_data.pressure_ref);

fprintf(fid, 'DIMENSIONLESS NUMBERS:\n');
fprintf(fid, '  Thermal diffusivity (alpha): %.2e m^2/s\n', problem_data.thermal_diffusivity);
fprintf(fid, '  Rayleigh number (Ra): %.2e\n', problem_data.rayleigh_number);
fprintf(fid, '    Formula: Ra = (rho_0 * beta * g * K * Delta_T * H) / (mu * alpha)\n');
fprintf(fid, '    where alpha = kappa_eff / (rho_0 * c_f)\n');
if problem_data.rayleigh_number < 40
    fprintf(fid, '    Flow regime: Conduction-dominated (Ra < 40)\n');
    fprintf(fid, '    Expected behavior: Minimal convection, heat transfer by conduction\n');
elseif problem_data.rayleigh_number < 100
    fprintf(fid, '    Flow regime: Transitional (40 < Ra < 100)\n');
    fprintf(fid, '    Expected behavior: Onset of convective instability\n');
else
    fprintf(fid, '    Flow regime: Convection-dominated (Ra > 100)\n');
    fprintf(fid, '    Expected behavior: Strong convective circulation\n');
end
fprintf(fid, '    Note: Critical Ra for porous media convection ~ 40 (Horton-Rogers-Lapwood)\n\n');

fprintf(fid, 'INITIAL CONDITIONS:\n');
fprintf(fid, '  Perturbation amplitude: %.2e K\n', problem_data.perturbation_amplitude);
fprintf(fid, '  Perturbation modes: %d\n\n', problem_data.perturbation_modes);

fprintf(fid, 'TIME INTEGRATION:\n');
fprintf(fid, '  Time step size: %.2f years (%.2e s)\n', problem_data.time_step / year, problem_data.time_step);
fprintf(fid, '  Final time: %.2f years (%.2e s)\n', problem_data.final_time / year, problem_data.final_time);
fprintf(fid, '  Expected time steps: %d\n', n_time_steps);
fprintf(fid, '  Actual time steps completed: %d\n', step);
fprintf(fid, '  Final simulation time: %.2f years (%.2e s)\n\n', time.t / year, time.t);

fprintf(fid, 'SOLVER SETTINGS:\n');
fprintf(fid, '  Newton max iterations: %d\n', problem_data.newton_max_iter);
fprintf(fid, '  Newton tolerance (update): %.2e\n', problem_data.newton_tol_du);
fprintf(fid, '  Newton tolerance (residual): %.2e\n', problem_data.newton_tol_res);
fprintf(fid, '  Print interval: %d steps\n\n', problem_data.print_interval);

fprintf(fid, 'FINAL SOLUTION STATISTICS:\n');
fprintf(fid, '  Temperature range: [%.2f, %.2f] K\n', min(T_final), max(T_final));
fprintf(fid, '  Mean temperature: %.2f K\n', mean(T_final));
fprintf(fid, '  Temperature std dev: %.2f K\n', std(T_final));
fprintf(fid, '  Pressure range: [%.2e, %.2e] Pa\n', min(p_final), max(p_final));
fprintf(fid, '  Mean pressure: %.2e Pa\n', mean(p_final));
fprintf(fid, '  Max velocity (x): %.2e m/s\n', max(abs(u_x)));
fprintf(fid, '  Max velocity (z): %.2e m/s\n', max(abs(u_z)));
fprintf(fid, '  Max |velocity|: %.2e m/s\n\n', max(sqrt(u_x.^2 + u_z.^2)));

fprintf(fid, 'OUTPUT FILES:\n');
fprintf(fid, '  - pressure_evolution_t_XXXXX.png (%d snapshots)\n', n_snapshots);
fprintf(fid, '  - temperature_evolution_t_XXXXX.png (%d snapshots)\n', n_snapshots);
fprintf(fid, '  - deltaT_evolution_t_XXXXX.png (%d snapshots)\n', n_snapshots);
fprintf(fid, '  - velocity_field_t_XXXXX.png (%d snapshots)\n', n_snapshots);
fprintf(fid, '  - problem_data.mat\n');
fprintf(fid, '  - simulation_report.txt\n');
fprintf(fid, '=======================================================\n');

fclose(fid);
fprintf('Saved: simulation_report.txt\n');

fprintf('\nVisualization complete.\n');
fprintf('All results saved to: %s\n', output_folder);

%% =======================================================================
%  HELPER FUNCTIONS
%% =======================================================================

function [residual, Jacobian] = applyBC(residual, Jacobian, y, grid, bc)
    % Apply Dirichlet boundary conditions to residual and Jacobian

    N = grid.Nx * grid.Nz;

    % Bottom temperature BC (T = T_hot)
    T_idx_bottom = 2 * bc.bottom_cells;
    residual(T_idx_bottom) = y(T_idx_bottom) - bc.T_hot;
    Jacobian(T_idx_bottom, :) = 0;
    Jacobian(T_idx_bottom, T_idx_bottom) = speye(length(T_idx_bottom));

    % Top temperature BC (T = T_cold)
    T_idx_top = 2 * bc.top_cells;
    residual(T_idx_top) = y(T_idx_top) - bc.T_cold;
    Jacobian(T_idx_top, :) = 0;
    Jacobian(T_idx_top, T_idx_top) = speye(length(T_idx_top));

    % Reference pressure BC
    p_idx_ref = 2 * bc.p_ref_cell - 1;
    residual(p_idx_ref) = y(p_idx_ref) - bc.p_ref;
    Jacobian(p_idx_ref, :) = 0;
    Jacobian(p_idx_ref, p_idx_ref) = 1;
end

function [u_x, u_z] = computeVelocity(y, grid, rock, fluid, physics, ...
                                       connectionsx, connectionsz)
    % Compute cell-centered velocity from pressure and temperature

    N = grid.Nx * grid.Nz;

    p = y(1:2:end);
    T = y(2:2:end);

    % Compute density
    rho = fluid.rho_0 * (1 - fluid.beta * (T - physics.T0));

    % Initialize velocities
    u_x = zeros(N, 1);
    u_z = zeros(N, 1);

    % X-direction velocities
    for f = 1:size(connectionsx, 1)
        i = connectionsx(f, 1);
        j = connectionsx(f, 2);

        % Harmonic mean permeability
        k_face = 2 / (1/rock.kx(i) + 1/rock.kx(j));

        % Pressure gradient
        dp_dx = (p(j) - p(i)) / grid.dx;

        % Darcy velocity
        u_face = -(k_face / fluid.mu) * dp_dx;

        % Distribute to cells
        u_x(i) = u_x(i) + u_face * 0.5;
        u_x(j) = u_x(j) + u_face * 0.5;
    end

    % Z-direction velocities
    for f = 1:size(connectionsz, 1)
        i = connectionsz(f, 1);
        j = connectionsz(f, 2);

        % Harmonic mean permeability
        k_face = 2 / (1/rock.kz(i) + 1/rock.kz(j));

        % Average density at face
        rho_face = 0.5 * (rho(i) + rho(j));

        % Pressure gradient with gravity
        dp_dz = (p(j) - p(i)) / grid.dz;
        gravity_term = rho_face * physics.g;

        % Darcy velocity (z is positive upward)
        u_face = -(k_face / fluid.mu) * (dp_dz - gravity_term);

        % Distribute to cells
        u_z(i) = u_z(i) + u_face * 0.5;
        u_z(j) = u_z(j) + u_face * 0.5;
    end
end

