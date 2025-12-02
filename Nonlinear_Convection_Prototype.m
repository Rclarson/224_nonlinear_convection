% Energy 224: Advanced Subsurface Flow Simulations
% Final Project
% Richard Larson

% Part A: Intialization

% Inputs

clc; clear; close all;


grid.p_open_bc = false;
grid.T_open_bc = false;

grid.Lx = 1;
grid.Ly = 1;
grid.Lz = 1;
% grid.Dtop = 5000;
% grid.Dbot = 5000;

grid.Nx = 30;
grid.Ny = 1;
grid.Nz = 30;


test_variable = 1;
rock.K_heat = ones(grid.Nx*grid.Nz);% * test_variable;
Ke = .2* test_variable ;  % should be 0.6 in correct units

rock.kx = ones(grid.Nx*grid.Nz,1 ) * test_variable*1;
rock.phi = ones(grid.Nx*grid.Nz, 1) ;
rock.kz = ones(grid.Nx*grid.Nz, 1) * test_variable;

T_init_bc = 1;
fluid.mu = 1;



% Calculated initialziations

grid.dx = grid.Lx / grid.Nx;
grid.dy = grid.Ly / grid.Ny;
grid.dz = grid.Lz / grid.Nz;


connectionsx = connectionsListx(grid);
connectionsz = connectionsListz(grid);


% boundary condition initialzations
xl = 1;
yl = 1;
xyl = 1;


connections_bc_x=[];
connections_bc_z=[];




if grid.p_open_bc == false
    
elseif grid.p_open_bc == true
    rock.kx(end+1) = rock.kx(end);
    rock.kz(end+1) = rock.kz(end);
    rock.phi(end+1) = rock.phi(end);
end


if grid.T_open_bc == false
    T_init_bc = ones(grid.Nx*grid.Nz,1) * T_init_bc;
elseif grid.T_open_bc == true
    T_init_bc = ones(grid.Nx*grid.Nz + 1,1) * T_init_bc;
end



t_final = 300;
grid.t_end = t_final;
dt = 10;
grid.dt = dt;
grid.max_dt = 1;
Nt = t_final/dt;

Tn = T_init_bc;
step = 1;
T_empty = 1;
t_step = 1;
Tn_ot = [];
T0 = Tn;
T = Tn;
rho_0 = 1;
alpha = 1;
c  = 1;

[Kex, Kez] = computeConductivity(grid, rock, fluid, connectionsx, connectionsz, T_empty, Ke);

rho = boussinesqDensity(T, T0, rho_0, alpha);

if grid.p_open_bc == false
    [Tx, Tz] = transmissibility(grid, rock, fluid, connectionsx, connectionsz, rho) ;
elseif grid.p_open_bc == true
    [Tx, Tz] = transmissibility(grid, rock, fluid, connectionsx_tot, connectionsz_tot, rho) ;
end
    
%GammaA = sparse(computeGammaA(grid,rock,fluid, connectionsx, connectionsz,connections_bc_x,connections_bc_z, Tx, Ty, rho));

%H = computeEnthalpy(grid, rock,dt, T,T0, c);


t=0;
t_initial = 0;
yv1 = ones(grid.Nx*grid.Ny*grid.Nz*2,1);
y = ones(grid.Nx*grid.Ny*grid.Nz*2,1);
s=1;

e1_crit = 1;
e2_crit = 1000;
e3_crit = 1;

while t < grid.t_end
    t_steps(s) =  t;
    y_solution(:,s) = yv1;
    
    s = s + 1;
    if t == t_initial
        yn = y;
        % yn(end+1)=0
        % yn(end+1)=0
    else
        yn = yv1;
        %ynp1 = ynp1;
    end
    epsilon = 1;
    e1 = 1;
    e2 = 1;
    e3 = 1;
    iter = 0;

    convergence = 0;
        %while epsilon >= epsilon_input
        while convergence == 0 && iter < 40 
            %rename connectionsz at some point
            residual = computeResidual(grid, rock, fluid, connectionsx, connectionsz, ...
                                    y, yv1, dt, T_empty, Ke, ...
                                    T0, rho_0, alpha, c);

            Jacobian = computeJacobian(grid, rock, fluid, connectionsx, connectionsz, ...
                        y, yv1, dt, T_empty, Ke, ...
                        T0, rho_0, alpha, c);

            % Temperature BC attempts
            N = grid.Nx*grid.Ny*grid.Nz;
            N_start = N - grid.Nx+1;
            bc_cells = [N_start:1:N];   % column vector
            T_idx    = 2 * bc_cells;         % column

            T_bc_val = 1.0;                    % scalar BC, same for all
            T_bc     = T_bc_val * ones(numel(T_idx),1);   % column, length = numel(T_idx)

            % Residual BC
            residual(T_idx) = yv1(T_idx) - T_bc;

            % Jacobian BC
            Jacobian(T_idx, :)     = 0;
            Jacobian(T_idx, T_idx) = speye(numel(T_idx));







            N = grid.Nx*grid.Ny*grid.Nz;
            N_start = N - grid.Nx+1;
            bc_cells = [445];   % column vector
            T_idx    = 2 * bc_cells;         % column

            T_bc_val = 1-.0000001;                    % scalar BC, same for all
            T_bc     = T_bc_val * ones(numel(T_idx),1);   % column, length = numel(T_idx)

            % Residual BC
            residual(T_idx) = yv1(T_idx) - T_bc;

            % Jacobian BC
            Jacobian(T_idx, :)     = 0;
            Jacobian(T_idx, T_idx) = speye(numel(T_idx));







    

            % pressure BC attempts
            N = grid.Nx*grid.Ny*grid.Nz;
            N_start = N - grid.Nx+1;
            bc_cells = [N_start:1:N];   % column vector
            P_idx    = 2 * bc_cells - 1;         % column

            P_bc_val = 1;                    % scalar BC, same for all
            P_bc     = P_bc_val * ones(numel(P_idx),1);   % column, length = numel(T_idx)

            % Residual BC
            residual(P_idx) = yv1(P_idx) - P_bc;

            % Jacobian BC
            Jacobian(P_idx, :)     = 0;
            Jacobian(P_idx, P_idx) = speye(numel(P_idx));
            % % 
            % 
            % 
            % %---  
            % p_ref_cell = 21;
            % p_ref_val  = 1;
            % p_idx      = 2 * p_ref_cell - 1;
            % 
            % residual(p_idx) = yv1(p_idx) - p_ref_val;
            % Jacobian(p_idx, :)    = 0;
            % Jacobian(p_idx, p_idx) = 1;
            % 
            % 
            % Q = zeros(grid.Nx*grid.Ny*grid.Nz*2, 1);
            % 
            % Q((1:1:grid.Nx)*2) = 0.0000001;
            % 
            % %Q((1:1:grid.Nx)*2-1) = -0.000000001;
            % Q((110)*2) = 0.0000001;
            % 
            % if s <= 2
            % 
            %     residual = residual - Q;
            % 
            % end





            dy = -Jacobian \ residual;   % Newton update
        
            %ynp1 = yn + dy;        
            yv1 = yv1 + dy;   
            epsilon = norm(dy);  

            iter = iter + 1;

            e1 = 1;
            
            %e2 = max(abs(dy(2:2:end)))
            e2 = norm(dy)

            
            e3 = 1;
            


            if e1 <= e1_crit && e2 <= e2_crit && e3 <= e3_crit
                convergence = 1;
            end
        end 




        t = t + grid.dt;
        
        

end






pmap = y_solution(1:2:end,end);
Tmap = y_solution(2:2:end,end);
% 
% imagesc( flipud(transpose(reshape(pmap, grid.Nx, grid.Ny))))
% set(gca, 'FontSize', 14);
% title('density');
% xlabel('x');
% ylabel('y');
% c = colorbar;
% set(c, 'FontSize', 14);
% colormap jet

figure 
imagesc( flipud(transpose(reshape(Tmap, grid.Nx, grid.Nz))))
set(gca, 'FontSize', 14);
title('Temperature');
xlabel('x');
ylabel('y');
c = colorbar;
set(c, 'FontSize', 14);
colormap jet


