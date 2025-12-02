function rho = boussinesqDensity(T, T0, rho_0, alpha)
    T0=T0(1);
    T0 = T0.*ones(length(T),1);
    rho = rho_0*(1-alpha*(T-T0));
end