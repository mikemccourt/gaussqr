function fobj = ex15_gqr_resnorm(coef,GQR,x,uold,dt,BC,t)
    usol = @(x,t) exp(-t)*(1-x.^2);
    usolx = @(x,t) exp(-t)*(-2*x);
    
    GQR.coef = coef; % Substitute in the test coefficients
    fu = ex15_gqr_residual(GQR,x,uold,dt);
    
    BCenforce = 1/dt;
%     BCenforce = 1;
    
    switch BC(1)
        case 0 % Dirichlet
            fu(1) = gqr_eval(GQR,x(1))-usol(x(1),t);
        case 1 % Neumann
            fu(1) = gqr_eval(GQR,x(1),1)-usolx(x(1),t);
    end
    
    switch BC(end)
        case 0 % Dirichlet
            fu(end) = gqr_eval(GQR,x(end))-usol(x(end),t);
        case 1 % Neumann
            fu(end) = gqr_eval(GQR,x(end),1)-usolx(x(end),t);
    end
    
    fu([1,end]) = fu([1,end])*BCenforce;
    
    fobj = norm(fu);
    fobj = fu;
end