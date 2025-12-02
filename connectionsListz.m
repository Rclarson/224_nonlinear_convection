function connectionsz = connectionsListz(grid)
    s = 1;
    for j = 1:grid.Nz-1
        for i = 1:grid.Nx
            l = (j-1)* grid.Nx +i;
            connectionsz(s,1)= l;
            connectionsz(s,2)= l + grid.Nz;
            s = s+1;
        end
    end
end