function [Tx, Tz] = transmissibility(grid, rock, fluid, connectionsx, connectionsz,rho)       
        connection_kx = 2./(1./rock.kx(connectionsx(:,1))+1./rock.kx(connectionsx(:,2)));

        % for upstreaming can I just take the larger value of rho?  but
        % then it doesn't work so well with horizontal and potentials
        Bavg = (rho(connectionsx(:,1)) + rho((connectionsx(:,2))))./2;
        Bavg = 1;
        T = (1./(Bavg.*fluid.mu).*(grid.dy.*grid.dz.*connection_kx./grid.dx)); % 887.5 or not
        Tx = T;
        connection_kz = 2./(1./rock.kz(connectionsz(:,1))+1./rock.kz(connectionsz(:,2)));
    
        Bavg = (rho(connectionsz(:,1)) + rho((connectionsz(:,2))))./2;
        Bavg = 1;
        T = (1./(Bavg.*fluid.mu).*(grid.dx.*grid.dy.*connection_kz./grid.dz)); % 887.5 or not
        Tz = T;
end