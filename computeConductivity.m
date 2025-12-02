function [Kex, Kez] = computeConductivity(grid, rock, fluid, connectionsx, connectionsz, T, Ke)       
        Kex = Ke*ones(length(connectionsx),1)*(grid.dy.*grid.dz./grid.dx); % double check on this 
        Kez = Ke*ones(length(connectionsz),1)*(grid.dx.*grid.dy./grid.dz);
end

