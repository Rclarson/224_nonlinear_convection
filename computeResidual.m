function residual = computeResidualRobot(grid, rock, fluid, ...
                                          connectionsx, connectionsz, ...
                                          y, yv1, dt, T_empty, Ke, ...
                                          T0, rho_0, alpha, c)
% Residual for coupled p–T system with upwind ρ,h and gravity
%
% State vector ordering: [p1; T1; p2; T2; ...; pN; TN]

    %------------------------------------------------------------------
    % Basic sizes and indices
    %------------------------------------------------------------------
    N        = grid.Nx * grid.Ny * grid.Nz;
    p_idx    = 1:2:(2*N);     % pressure DOFs
    T_idx    = 2:2:(2*N);     % temperature DOFs

    % Make sure T0 is scalar (you’ve been using a uniform reference)
    T0 = T0(1);

    %------------------------------------------------------------------
    % Cell-wise primary variables and properties at t_{n+1}
    %------------------------------------------------------------------
    p_cells = yv1(p_idx);               % N x 1
    T_cells = yv1(T_idx);               % N x 1

    rho = boussinesqDensity(T_cells, T0, rho_0, alpha);    % N x 1
    h   = (T_cells - T0) .* c;                             % N x 1

    % Transmissibilities and effective conductivities (face-based)
    [Tx, Tz]   = transmissibility(grid, rock, fluid, ...
                                  connectionsx, connectionsz, rho);
    [Kex, Kez] = computeConductivity(grid, rock, fluid, ...
                                     connectionsx, connectionsz, ...
                                     T_empty, Ke);

    %------------------------------------------------------------------
    % Upwind selection on faces (based on pressure)
    %------------------------------------------------------------------
    % X-direction faces
    ix  = connectionsx(:,1);      % left cell index
    jx  = connectionsx(:,2);      % right cell index

    p1x = p_cells(ix);
    p2x = p_cells(jx);

    h1x   = h(ix);
    h2x   = h(jx);
    rho1x = rho(ix);
    rho2x = rho(jx);

    mask_ix = (p1x >= p2x);       % upwind = i where pi >= pj
    mask_jx = ~mask_ix;

    rho_face_x = rho1x .* mask_ix + rho2x .* mask_jx;   % upwind ρ
    h_face_x   = h1x   .* mask_ix + h2x   .* mask_jx;   % upwind h

    % Z-direction faces
    iz  = connectionsz(:,1);
    jz  = connectionsz(:,2);

    p1z = p_cells(iz);
    p2z = p_cells(jz);

    h1z   = h(iz);
    h2z   = h(jz);
    rho1z = rho(iz);
    rho2z = rho(jz);

    mask_iz = (p1z >= p2z);
    mask_jz = ~mask_iz;

    rho_face_z = rho1z .* mask_iz + rho2z .* mask_jz;
    h_face_z   = h1z   .* mask_iz + h2z   .* mask_jz;

    %------------------------------------------------------------------
    % Assemble connection matrix R in the usual (l1,l2,r1,r2) style
    %------------------------------------------------------------------
    R = sparse(2*N, 2*N);

    %% ------------ X–direction blocks -------------------------------
    l1 = connectionsx(:,1)*2 - 1;   % p_i rows
    l2 = connectionsx(:,1)*2;       % T_i rows
    r1 = connectionsx(:,2)*2 - 1;   % p_j rows
    r2 = connectionsx(:,2)*2;       % T_j rows

    l1in = sub2ind(size(R), l1, r1);
    l2in = sub2ind(size(R), l2, r1);
    l3in = sub2ind(size(R), l1, r2);
    l4in = sub2ind(size(R), l2, r2);

    r1in = sub2ind(size(R), r1, l1);
    r2in = sub2ind(size(R), r2, l1);
    r3in = sub2ind(size(R), r1, l2);
    r4in = sub2ind(size(R), r2, l2);

    ld1in = sub2ind(size(R), l1, l1);
    ld2in = sub2ind(size(R), l2, l1);
    ld3in = sub2ind(size(R), l1, l2);
    ld4in = sub2ind(size(R), l2, l2);

    rd1in = sub2ind(size(R), r1, r1);
    rd2in = sub2ind(size(R), r2, r1);
    rd3in = sub2ind(size(R), r1, r2);
    rd4in = sub2ind(size(R), r2, r2);

    % Mass equation (p-rows) – advective flux with upwind ρ
    R(l1in)  = Tx .* rho_face_x;     % coeff for p_j in eqn of p_i
    R(r1in)  = Tx .* rho_face_x;     % symmetric
    R(ld1in) = R(ld1in) - Tx .* rho_face_x;
    R(rd1in) = R(rd1in) - Tx .* rho_face_x;

    % Energy equation (T-rows) – advective + conductive in x
    R(l2in)  = Tx .* rho_face_x .* h_face_x;
    R(r2in)  = Tx .* rho_face_x .* h_face_x;
    R(ld2in) = R(ld2in) - Tx .* rho_face_x .* h_face_x;
    R(rd2in) = R(rd2in) - Tx .* rho_face_x .* h_face_x;

    % Conduction part
    R(l4in)  = Kex;
    R(r4in)  = Kex;
    R(ld4in) = R(ld4in) - Kex;
    R(rd4in) = R(rd4in) - Kex;

    %% ------------ Z–direction blocks -------------------------------
    l1 = connectionsz(:,1)*2 - 1;
    l2 = connectionsz(:,1)*2;
    r1 = connectionsz(:,2)*2 - 1;
    r2 = connectionsz(:,2)*2;

    l1in = sub2ind(size(R), l1, r1);
    l2in = sub2ind(size(R), l2, r1);
    l3in = sub2ind(size(R), l1, r2);
    l4in = sub2ind(size(R), l2, r2);

    r1in = sub2ind(size(R), r1, l1);
    r2in = sub2ind(size(R), r2, l1);
    r3in = sub2ind(size(R), r1, l2);
    r4in = sub2ind(size(R), r2, l2);

    ld1in = sub2ind(size(R), l1, l1);
    ld2in = sub2ind(size(R), l2, l1);
    ld3in = sub2ind(size(R), l1, l2);
    ld4in = sub2ind(size(R), l2, l2);

    rd1in = sub2ind(size(R), r1, r1);
    rd2in = sub2ind(size(R), r2, r1);
    rd3in = sub2ind(size(R), r1, r2);
    rd4in = sub2ind(size(R), r2, r2);

    % Mass equation (p-rows) – vertical advective flux (no gravity here;
    % gravity goes into G below)
    R(l1in)  = Tz .* rho_face_z;
    R(r1in)  = Tz .* rho_face_z;
    R(ld1in) = R(ld1in) - Tz .* rho_face_z;
    R(rd1in) = R(rd1in) - Tz .* rho_face_z;

    % Energy equation (T-rows) – vertical advective + conductive should
    % this also have c or alpha?
    R(l2in)  = Tz .* rho_face_z .* h_face_z;
    R(r2in)  = Tz .* rho_face_z .* h_face_z;
    R(ld2in) = R(ld2in) - Tz .* rho_face_z .* h_face_z;
    R(rd2in) = R(rd2in) - Tz .* rho_face_z .* h_face_z;
    % 
    % gamma      = computeGamma(rho, fluid);   % N x 1
    % gamma_avg  = -0.5 * (gamma(iz) + gamma(jz));
    % 
    % R(l3in)  = Tz .* gamma_avg;
    % R(r3in)  = -Tz .* gamma_avg;
    % R(ld3in) = R(ld3in) - Tz .* gamma_avg;
    % R(rd3in) = R(rd3in) + Tz .* gamma_avg;

    R(l4in)  = Kez;
    R(r4in)  = Kez;
    R(ld4in) = R(ld4in) - Kez;
    R(rd4in) = R(rd4in) - Kez;

    %------------------------------------------------------------------
    % Gravity source term G (acts on p- and T-equations in z)
    %------------------------------------------------------------------
    G = sparse(2*N, 1);

    gamma      = computeGamma(rho, fluid);   % N x 1
    gamma_avg  = 0.5 * (gamma(iz) + gamma(jz));

    % NOTE: this uses the same "overwrite / last-write wins" pattern
    % you’ve used before; for a cleaner sum you’d use accumarray.
    G(2*iz-1) = G(2*iz-1) - Tz .* gamma_avg .* grid.dz;
    G(2*jz-1) = G(2*jz-1) + Tz .* gamma_avg .* grid.dz;

    G(2*iz)   = G(2*iz)   - h_face_z .* Tz .* gamma_avg .* grid.dz;
    G(2*jz)   = G(2*jz)   + h_face_z .* Tz .* gamma_avg .* grid.dz;

    %------------------------------------------------------------------
    % Accumulation term B (energy only, no compressibility in mass)
    %------------------------------------------------------------------
    B = sparse(2*N, 1);

    T_prev = y(T_idx);         % T^n
    T_new  = yv1(T_idx);       % T^{n+1}

    accum_n_p  = (grid.dx*grid.dy*grid.dz/dt) .* boussinesqDensity(T_prev,T0,rho_0,alpha);
    accum_np_p = (grid.dx*grid.dy*grid.dz/dt) .* boussinesqDensity(T_prev,T0,rho_0,alpha);

    accum_n  = computeEnthalpy(grid, rock, dt, T_prev, T0, c).* boussinesqDensity(T_prev,T0,rho_0,alpha);
    accum_np = computeEnthalpy(grid, rock, dt, T_new,  T0, c).* boussinesqDensity(T_new,T0,rho_0,alpha);

    % Only T-rows get accumulation
    %B(p_idx) = accum_np_p - accum_n_p;
    B(T_idx) = accum_np - accum_n;

    %------------------------------------------------------------------
    % Final residual
    %------------------------------------------------------------------
    residual = R * yv1 + B - G;
end
