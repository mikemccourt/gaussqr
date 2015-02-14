function [ret1,ret2] = SurrModel4_CV(par,rbf,x,y,cvind)
% This function is required for the fmincon solve to optimize ep and mu
% It can be called in two ways: to evaluate the CV residual, or constraint
% Don't try and call it with the nonlinear constraint though
% Eventually will be called with rbfs = {rbf,rbfdx}
%
% The way this should be called is
%   opt_par = fmincon(@(epmu)SurrModel4_CV(epmu,rbfs,x,y,cvind),epmu_guess,[],[],[],[],[0;0],[Inf;Inf])
% epmu and epmu_guess should be a column vector [ep;mu]
% Or, if you use the unconstrained version, you need 
%   opt_par = fminunc(@(epmu)SurrModel4_CV(exp(epmu),rbfs,x,y,cvind),log(epmu_guess))
% This is basically a rescaling of the problem to the log domain
%
% CVres = SurrModel4_CV(par,rbf,x,y,cvind)
%   Inputs:
%      par - [ep;mu]
%      rbf - kernel for integration: rbf(DistanceMatrix(x,z,epvec))
%      x - vector of N data location
%      y - vector of N data values at locations x
%      cvind - {out_1,...,out_k} with k sets of indices left out
%   Output:
%      CVres - the CV residual, summed over the cvind
%              the norm type is given in GAUSSQR_PARAMETERS
%
%%%%%
% The constraint stuff is not really working right now
% I'm going to figure it out eventually
% 
% [c,ceq] = SurrModel4_CV(par,rbfs,x,y,xtest)
%   Outputs:
%      c - The inequality constraint, computed as the sum of
%             -sign(k(x)^T(K+mu*I)y), at xtest for CDF>0
%             
%               
%      ceq - The equality constraint, not used for this method
%                this will always return 0

ep = par(1);
mu = par(2);
N = size(x,1);
switch nargout
    case {0,1}
        Ncv = length(cvind);
        
        residual = 0;
        for k=1:Ncv
            this_out = cvind{k};
            this_in = setdiff(1:N,this_out);
            xout = x(this_out);
            yout = y(this_out);
            xin = x(this_in);
            yin = y(this_in);
            Kin = rbf(DistanceMatrix(xin,xin,ep));
            I = eye(size(xin,1));
            Kout = rbf(DistanceMatrix(xout,xin,ep));
            this_res = norm(Kout*((Kin + mu*I)\yin)-yout);
            residual = residual + this_res^2;
        end
        ret1 = residual;
    case 2
        error('This constraint business has not yet been implemented')
    otherwise
        error('You can''t get more than two outputs from this function')
end






