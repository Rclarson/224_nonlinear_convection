function [Tx, Tz] = transmissibility(grid, rock, fluid, connectionsx, connectionsz, rho)
    % X-direction transmissibilities
    connection_kx = 2./(1./rock.kx(connectionsx(:,1)) + 1./rock.kx(connectionsx(:,2)));

    % For Boussinesq, we typically use constant density in transmissibility
    % (density variation only in buoyancy term)
    % Using rho_0 = 1 in dimensionless form, or can use actual rho_0
    Bavg_x = ones(size(connection_kx));  % Constant for Boussinesq
    Tx = (1./(Bavg_x.*fluid.mu) .* (grid.dy.*grid.dz.*connection_kx./grid.dx));

    % Z-direction transmissibilities
    connection_kz = 2./(1./rock.kz(connectionsz(:,1)) + 1./rock.kz(connectionsz(:,2)));

    Bavg_z = ones(size(connection_kz));  % Constant for Boussinesq
    Tz = (1./(Bavg_z.*fluid.mu) .* (grid.dx.*grid.dy.*connection_kz./grid.dz));
end