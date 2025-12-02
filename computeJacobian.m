function Jacobian = computeJacobian(grid, rock, fluid, connectionsx, connectionsz, ...
                                    y, yv1, dt, T_empty, Ke, ...
                                    T0, rho_0, alpha, c)
    
    N = grid.Nx*grid.Ny*grid.Nz;
    Jacobian = sparse(N*2, N*2);

    m = 1:2:N*2;
    l = 2:2:N*2;

    % do we really need upwinding for these equations
    % also need to check the yv1-yv1 for sign, might be different in second
    % block
    rho = boussinesqDensity(yv1(l),T0,rho_0,alpha);
    

    [Tx, Tz] = transmissibility(grid, rock, fluid, connectionsx, connectionsz, rho);

    [Kex, Kez] = computeConductivity(grid, rock, fluid, connectionsx, connectionsz, T_empty, Ke);

    drhodt = computeRhoDerivative(T0(1),rho_0,alpha);
    dhdt = computeEnthalpyDerivative(T0(1),c);

    
    rho_avg_x = (boussinesqDensity(yv1(connectionsx(:,1)*2),T0,rho_0,alpha)+boussinesqDensity(yv1(connectionsx(:,2)*2),T0,rho_0,alpha))/2;
    rho_avg_z = (boussinesqDensity(yv1(connectionsz(:,1)*2),T0,rho_0,alpha)+boussinesqDensity(yv1(connectionsz(:,2)*2),T0,rho_0,alpha))/2;

    

    h = (yv1(l)-T0).*c;


    % using averages now but may need to change for correctness

    p = yv1(1:2:end);

    i = connectionsx(:,1);   % left cell
    j = connectionsx(:,2);   % right cell
    
    p1 = p(i);
    p2 = p(j);
    
    h1 = h(i);
    h2 = h(j);

    rho1 = rho(i);
    rho2 = rho(j);
    
    % Upwind selection: enthalpy at the higher-pressure cell
    h_upstream_x = h1 .* (p1 >= p2) + h2 .* (p2 > p1);
    rho_upstream_x = rho1 .* (p1 >= p2) + rho2 .* (p2 > p1);

    i = connectionsz(:,1);   % left cell
    j = connectionsz(:,2);   % right cell
    
    p1 = p(i);
    p2 = p(j);
    
    h1 = h(i);
    h2 = h(j);

    rho1 = rho(i);
    rho2 = rho(j);
        
    % Upwind selection: enthalpy at the higher-pressure cell
    h_upstream_z = h1 .* (p1 >= p2) + h2 .* (p2 > p1);
    rho_upstream_z = rho1 .* (p1 >= p2) + rho2 .* (p2 > p1);

    
    
    % h_upstream_x = max([p(connectionsx(:,1)) p(connectionsx(:,2))], 2 );
    % h_upstream_z = max([p(connectionsz(:,1)) p(connectionsz(:,2))], 2 );

    h_avg_x = (h(connectionsx(:,1)) + h(connectionsx(:,2)))/2;
    h_avg_z = (h(connectionsz(:,1)) + h(connectionsz(:,2)))/2;

    h_avg_x = h_upstream_x;
    h_avg_z = h_upstream_z;

    rho_avg_x = rho_upstream_x;
    rho_avg_z = rho_upstream_z;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p_cells = yv1(m);
    
    ix = connectionsx(:,1);   % left cell
    jx = connectionsx(:,2);   % right cell
    
    p1 = p_cells(ix);
    p2 = p_cells(jx);
    
    h1 = h(ix);   h2 = h(jx);
    rho1 = rho(ix); rho2 = rho(jx);
    
    mask_i = (p1 >= p2);   % upstream = i?
    mask_j = ~mask_i;
    
    rho_face = rho1 .* mask_i + rho2 .* mask_j;
    h_face   = h1   .* mask_i + h2   .* mask_j;
    
    drho_dT_i = drhodt .* mask_i;   % scalar drho_dT_j times mask
    drho_dT_j = drhodt.* mask_j;
    dh_dT_i   = dhdt   .* mask_i;
    dh_dT_j   = dhdt   .* mask_j;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5



    
    %if p1 >= p2


    l1 = connectionsx(:,1)*2 - 1;
    l2 = connectionsx(:,1)*2;
    r1 = connectionsx(:,2)*2 - 1;
    r2 = connectionsx(:,2)*2;

    l1in = sub2ind(size(Jacobian), l1, r1);  
    l2in = sub2ind(size(Jacobian), l2, r1);  
    l3in = sub2ind(size(Jacobian), l1, r2);  
    l4in = sub2ind(size(Jacobian), l2, r2);


    r1in = sub2ind(size(Jacobian), r1, l1);  
    r2in = sub2ind(size(Jacobian), r2, l1);  
    r3in = sub2ind(size(Jacobian), r1, l2);  
    r4in = sub2ind(size(Jacobian), r2, l2);  



    ld1in = sub2ind(size(Jacobian), l1, l1);  
    ld2in = sub2ind(size(Jacobian), l2, l1);  
    ld3in = sub2ind(size(Jacobian), l1, l2);  
    ld4in = sub2ind(size(Jacobian), l2, l2);  

    rd1in = sub2ind(size(Jacobian), r1, r1);  
    rd2in = sub2ind(size(Jacobian), r2, r1);  
    rd3in = sub2ind(size(Jacobian), r1, r2);  
    rd4in = sub2ind(size(Jacobian), r2, r2);      


    dp_x = -(yv1(connectionsx(:,1)*2-1)-yv1(connectionsx(:,2)*2-1));


    Jacobian(l1in) = Tx .* rho_avg_x;  % this should be good for now
    Jacobian(l3in) = Tx .* drho_dT_i .* (dp_x);  % make sure indexing correctly;  adding a term for the gamma contributions (zero in x)
    Jacobian(l2in) = Tx .* rho_avg_x .* h_avg_x;  % this should be good for now
    Jacobian(l4in) = Tx .* (dh_dT_i .*  rho_avg_x + drho_dT_i  .* h_avg_x) .* (dp_x)  - Kex; 

    Jacobian(ld1in) = Jacobian(ld1in) - Tx .* rho_avg_x;
    Jacobian(ld3in) = Jacobian(ld3in) - Tx .* drho_dT_i .* ((dp_x));
    Jacobian(ld2in) = Jacobian(ld2in) - Tx .* rho_avg_x .* h_avg_x;
    Jacobian(ld4in) = Jacobian(ld4in) - (Tx .* (dh_dT_i .*  rho_avg_x + drho_dT_i  .* h_avg_x) .* ((dp_x)) - Kex); 

    % changed this section from the upwinded version

    Jacobian(r1in) = Tx .* rho_avg_x ;
    Jacobian(r3in) = Tx .* drho_dT_j .* (-(dp_x));  % make sure indexing correctly;  adding a term for the gamma contributions (zero in x)
    Jacobian(r2in) = Tx .* rho_avg_x .* h_avg_x;
    Jacobian(r4in) = Tx .* (dh_dT_j .*  rho_avg_x + drho_dT_j  .* h_avg_x) .* (-(dp_x))  - Kex; % need to add other part of the derivative


    Jacobian(rd1in) = Jacobian(rd1in) - Tx .* rho_avg_x;
    Jacobian(rd3in) = Jacobian(rd3in) - Tx .* drho_dT_j .* (-(dp_x));
    Jacobian(rd2in) = Jacobian(rd2in) - Tx .* rho_avg_x .* h_avg_x;
    Jacobian(rd4in) = Jacobian(rd4in) - (Tx .* (dh_dT_j .*  rho_avg_x + drho_dT_j  .* h_avg_x) .* (-(dp_x)) - Kex); 

    %end


    % % needs upwinding

    %%%%%%%%%%%%%%%%%%%%%%%%
    l1 = connectionsz(:,1)*2 - 1;
    l2 = connectionsz(:,1)*2;
    r1 = connectionsz(:,2)*2 - 1;
    r2 = connectionsz(:,2)*2;

    l1in = sub2ind(size(Jacobian), l1, r1);  
    l2in = sub2ind(size(Jacobian), l2, r1);  
    l3in = sub2ind(size(Jacobian), l1, r2);  
    l4in = sub2ind(size(Jacobian), l2, r2);


    r1in = sub2ind(size(Jacobian), r1, l1);  
    r2in = sub2ind(size(Jacobian), r2, l1);  
    r3in = sub2ind(size(Jacobian), r1, l2);  
    r4in = sub2ind(size(Jacobian), r2, l2);  



    ld1in = sub2ind(size(Jacobian), l1, l1);  
    ld2in = sub2ind(size(Jacobian), l2, l1);  
    ld3in = sub2ind(size(Jacobian), l1, l2);  
    ld4in = sub2ind(size(Jacobian), l2, l2);  

    rd1in = sub2ind(size(Jacobian), r1, r1);  
    rd2in = sub2ind(size(Jacobian), r2, r1);  
    rd3in = sub2ind(size(Jacobian), r1, r2);  
    rd4in = sub2ind(size(Jacobian), r2, r2);     


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p_cells = yv1(m);
    
    iz = connectionsz(:,1);   % left cell
    jz = connectionsz(:,2);   % right cell
    
    p1 = p_cells(iz);
    p2 = p_cells(jz);
    
    h1 = h(iz);   h2 = h(jz);
    rho1 = rho(iz); rho2 = rho(jz);
    
    mask_i = (p1 >= p2);   % upstream = i?
    mask_j = ~mask_i;
    
    rho_face = rho1 .* mask_i + rho2 .* mask_j;
    h_face   = h1   .* mask_i + h2   .* mask_j;
    
    drho_dT_i = drhodt .* mask_i;   % scalar drho_dT_j times mask
    drho_dT_j = drhodt .* mask_j;
    dh_dT_i   = dhdt   .* mask_i;
    dh_dT_j   = dhdt   .* mask_j;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5



    



    rho = boussinesqDensity(yv1(l),T0,rho_0,alpha);
    gamma = computeGamma(rho, fluid);

    %[gamma_w2, gamma_g2] = computeGamma(fluid , p2);

    dz_up = -grid.dz;
    dz_down = grid.dz;

    %gamma_avg = (gamma(connectionsz(:,1)) + gamma(connectionsz(:,2)))/2;
    %gamma_avg = 1*ones(length(connectionsz),1);
    gamma_avg = computeGamma(rho_upstream_z,fluid);

    dp_z = -(yv1(connectionsz(:,1)*2-1)-yv1(connectionsz(:,2)*2-1));


    Jacobian(l1in) = Tz .* rho_avg_z;
    Jacobian(l3in) = Tz .* drho_dT_i .* ((dp_z) + gamma_avg .* dz_up) - Tz .* drho_dT_i .* dz_up;  % check the gamma derivative
    

    Jacobian(rd1in) = Jacobian(rd1in) - Tz .* rho_avg_z;
    Jacobian(rd3in) = Jacobian(rd3in) - (Tz .* drho_dT_j .* (-(dp_z) + gamma_avg * dz_down) - Tz .* drho_dT_j  .* dz_down); % may need a sign change on grid.dz
    
    Jacobian(r1in) = Tz .* rho_avg_z;
    Jacobian(r3in) = Tz .* drho_dT_j .* (-(dp_z) + gamma_avg .* dz_down) - Tz .* drho_dT_j  .* dz_down;

    Jacobian(ld1in) = Jacobian(ld1in) - Tz .* rho_avg_z;
    Jacobian(ld3in) = Jacobian(ld3in) - (Tz .* drho_dT_i .* (dp_z + gamma_avg .* dz_up) - Tz .* drho_dT_i .* dz_up);





    % 

    Jacobian(l2in) = Tz .* rho_avg_z .* h_avg_z;
    Jacobian(l4in) = Tz .* (dh_dT_i .*  rho_avg_z + drho_dT_i  .* h_avg_z) .* (dp_z + gamma_avg*dz_up) - h_avg_z .* Tz .* drho_dT_i .* ( dz_up) - Kez; 

    Jacobian(rd2in) = Jacobian(rd2in) - Tz .* rho_avg_z .* h_avg_z;
    Jacobian(rd4in) = Jacobian(rd4in) - (Tz .* (dh_dT_j .*  rho_avg_z + drho_dT_j  .* h_avg_z) .* (-(dp_z) + gamma_avg*dz_down) - h_avg_z .* Tz .* drho_dT_j .* ( dz_down) - Kez); 

    Jacobian(r2in) = Tz .* rho_avg_z .* h_avg_z;
    Jacobian(r4in) = Tz .* (dh_dT_j .*  rho_avg_z + drho_dT_j  .* h_avg_z) .* (-(dp_z) + gamma_avg*dz_down) - h_avg_z .* Tz .* drho_dT_j .* ( dz_down) - Kez;

    Jacobian(ld2in) = Jacobian(ld2in) - Tz .* rho_avg_z .* h_avg_z;
    Jacobian(ld4in) = Jacobian(ld4in) - (Tz .* (dh_dT_i .*  rho_avg_z + drho_dT_i  .* h_avg_z) .* (dp_z + gamma_avg*dz_up) - h_avg_z .* Tz .* drho_dT_i .* ( dz_up) - Kez); 


    %Jacobian(m, m) = Jacobian(m, m) + (grid.dx*grid.dy*grid.dz/dt).* (rock.phi.* drho_dT_j);

    Jacobian(l, l) = Jacobian(l, l) - (grid.dx*grid.dy*grid.dz/dt).* (rock.phi.* drhodt .* h + rock.phi.* rho.* dhdt);


end