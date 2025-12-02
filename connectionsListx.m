function connectionsx = connectionsListx(grid)
    s = 1;
    for j = 1:grid.Nz
        for i = 1:grid.Nx - 1
            l = (j-1)* grid.Nx +i;
            connectionsx(s,1)= l;
            connectionsx(s,2)= l + 1;
            s = s+1;
        end
    end
end