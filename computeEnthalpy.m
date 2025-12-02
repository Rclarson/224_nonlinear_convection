function H = computeEnthalpy(grid, rock, dt, T, T0, c)

    h = (T-T0).*c; % from engineering toolbox c ~= 4.5 kJ/kg/C or 4.184 according to c_p google
    if grid.T_open_bc == false
        H = ones(grid.Nx*grid.Nz,1);
        %Hijk = grid.dx*grid.dy*grid.dz*rock.phi(1:length(T)).*h/dt;  when
        %using this commented line needs a different sign on diag(H)
        Hijk = grid.dx*grid.dy*grid.dz*rock.phi/dt;
        H = h.*H .* Hijk;        
    elseif grid.T_open_bc == true
        H = ones(grid.Nx*grid.Nz+1,1);
        %Hijk = grid.dx*grid.dy*grid.dz*rock.phi.*h/dt;
        Hijk = grid.dx*grid.dy*grid.dz*rock.phi./dt;
        H = h.*H .* Hijk; 
    end

end