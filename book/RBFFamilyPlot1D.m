% RBFFamilyPlot1D
%
% Generates plot of a family of RBF basic functions centered at origin
%
% Calls on: DistanceMatrix
%
% no error checking

function RBFFamilyPlot1D

%%%%%%%%%%%%%%%%%% Editable Section %%%%%%%%%%%%%%%%%%%%
% 
% Define the RBFs
%rbf1 = inline('max(1-epsilon*r,0).^2');    % Wendland C0
%rbf2 = inline('max(1-epsilon*r,0).^4.*(4*epsilon*r+1)');    % Wendland C2
%rbf3 = inline('max(1-epsilon*r,0).^6.*(35*(epsilon*r).^2+18*epsilon*r+3)/3');    % Wendland C4
%rbf4 = inline('max(1-epsilon*r,0).^8.*(32*(epsilon*r).^3+25*(epsilon*r).^2+8*epsilon*r+1)');    % Wendland C6
%rbf1 = 'Wendland2onehalf';    % missing Wendland (2,0.5)
%rbf2 = 'Wendland3threehalf';    % missing Wendland (3,1.5)
%rbf1 = inline('max(1-epsilon*r,0).^4.*(5*(epsilon*r).^3+20*(epsilon*r).^2+29*epsilon*r+16)/16');    % Wu(3,3) C0, SPD(R7) 
%rbf2 = inline('max(1-epsilon*r,0).^5.*(5*(epsilon*r).^4+25*(epsilon*r).^3+48*(epsilon*r).^2+40*epsilon*r+8)/8');    % Wu(2,3) C2, SPD(R5) 
%rbf3 = inline('max(1-epsilon*r,0).^6.*(5*(epsilon*r).^5+30*(epsilon*r).^4+72*(epsilon*r).^3+82*(epsilon*r).^2+36*epsilon*r+6)/6');    % Wu(1,3) C4, SPD(R3) 
%rbf4 = inline('max(1-epsilon*r,0).^7.*(5*(epsilon*r).^6+35*(epsilon*r).^5+101*(epsilon*r).^4+147*(epsilon*r).^3+101*(epsilon*r).^2+35*epsilon*r+5)/5');    % Wu(0,3) C6, SPD(R1) 
%rbf1 = 'buhmann';   % Buhmann C2
%rbf1 = inline('max(1-epsilon*r,0).^(7/2).*(1+7/2*epsilon*r-135/8*(epsilon*r).^2)');    % Gneiting (2,7)
%rbf2 = inline('max(1-epsilon*r,0).^5.*(1+5*epsilon*r-27*(epsilon*r).^2)');    % Gneiting (2,10)
%rbf3 = inline('max(1-epsilon*r,0).^(15/2).*(1+15/2*epsilon*r-391/8*(epsilon*r).^2)');    % Gneiting (2,15)
%rbf4 = inline('max(1-epsilon*r,0).^12.*(1+12*epsilon*r-104*(epsilon*r).^2)');    % Gneiting (2,24)
%rbf1 = inline('max(1-epsilon*r,0).^4.*(1+4*epsilon*r-15*(epsilon*r).^2)');    % Gneiting C2, R3
%rbf2 = inline('max(1-epsilon*r,0).^6.*(3+18*epsilon*r+3*(epsilon*r).^2-192*(epsilon*r).^3)/3');    % Gneiting C4, R3
%rbf3 = inline('max(1-epsilon*r,0).^8.*(15+120*epsilon*r+210*(epsilon*r).^2-840*(epsilon*r).^3-3465*(epsilon*r).^4)/15');    % Gneiting C6, R3
%rbf1 = inline('max(1-epsilon*r/2,0)');   % Euclid s=1
%rbf2 = inline('1/(2*pi)*(4*acos(epsilon*r/2)-epsilon*r.*sqrt(4-(epsilon*r).^2))');   % Euclid s=2
%rbf3 = inline('1-1/(32*pi)*((4 + 16*pi)*epsilon*r - (epsilon*r).^3)');   % Euclid s=3
%rbf4 = inline('2/pi*acos(epsilon*r/2)-(1/(32*pi))*sqrt(4-(epsilon*r).^2).*(20*epsilon*r-(epsilon*r).^3)');   % Euclid s=4
%rbf5 = inline('1 - 1/(64*pi^2)*((12+8*pi+32*pi^2)*epsilon*r-(3+2*pi)*(epsilon*r).^3)');   % Euclid s=5
%rbf1 = inline('8/15*exp(-(epsilon*r).^2).*(15/8-5/2*(epsilon*r).^2+0.5*(epsilon*r).^4)');     % Gaussian quadratic Laguerre
%rbf1 = inline('1./(1+(epsilon*r).^2)');   % generalized IMQ
%rbf2 = inline('(1/3)*(3-(epsilon*r).^2)./(1+(epsilon*r).^2).^3');   % linear generalized IMQ
%rbf3 = inline('(1/5)*(5-10*(epsilon*r).^2+(epsilon*r).^4)./(1+(epsilon*r).^2).^5');   % quadratic generalized IMQ
%rbf1 = inline('epsilon*r');    % Norm
%rbf1 = inline('exp(-epsilon*r)');    % Basic Matern
%rbf2 = inline('exp(-0.5*epsilon*r)');    % Basic Matern
%rbf3 = inline('exp(-0.25*epsilon*r)');    % Basic Matern
%rbf4 = inline('exp(-0.125*epsilon*r)');    % Basic Matern
%rbf5 = inline('exp(-0.6125*epsilon*r)');    % Basic Matern
%rbf1 = inline('exp(-epsilon*r) + epsilon*r');    % Tension spline
%rbf2 = inline('exp(-0.5*epsilon*r) + 0.5*epsilon*r');    % Tension spline
%rbf3 = inline('exp(-0.25*epsilon*r) + 0.25*epsilon*r');    % Tension spline
%rbf4 = inline('exp(-0.125*epsilon*r) + 0.125*epsilon*r');    % Tension spline
%rbf5 = inline('exp(-0.6125*epsilon*r) + 0.6125*epsilon*r');    % Tension spline
%rbf1 = inline('(epsilon*r).^3');    % Cubic
%rbf1 = inline('sqrt(2)*exp(-epsilon*r).*sin(epsilon*r + pi/4)');    % Sobolev spline (BTA)
%rbf2 = inline('exp(-0.5*epsilon*r).*sin(0.5*epsilon*r + pi/4)');    % Sobolev spline (BTA)
%rbf3 = inline('exp(-0.25*epsilon*r).*sin(0.25*epsilon*r + pi/4)');    % Sobolev spline (BTA)
%rbf4 = inline('exp(-0.125*epsilon*r).*sin(0.125*epsilon*r + pi/4)');    % Sobolev spline (BTA)
%rbf5 = inline('exp(-0.6125*epsilon*r).*sin(0.6125*epsilon*r + pi/4)');    % Sobolev spline (BTA)
%rbf1 = inline('cos(epsilon*r)'); % Cosine kernel
rbf1 = @(ep,r) exp(-ep*r).*(1+ep*r);    % C^2 Matern
rbf2 = @(ep,r) sqrt(2)*exp(-ep*r).*sin(ep*r + pi/4);    % H^2 Sobolev spline (BTA)

% Parameter for basis function
epsilon = 5;
%
% Resolution of evaluation grid for error computation and plotting
neval = 101;      % to create neval evaluation grid in [-1,1]
%
% Create neval equally spaced evaluation locations in [-1,1]
evalpoints = linspace(-2,2,neval)';
DM_eval = DistanceMatrix(evalpoints,0);
basisfunction1 = rbf1(epsilon,DM_eval);
basisfunction2 = rbf2(3,DM_eval);
%basisfunction3 = feval(rbf1,2*epsilon,DM_eval);
%basisfunction2 = feval(rbf2,epsilon,DM_eval);
%basisfunction3 = feval(rbf3,epsilon,DM_eval);
%basisfunction4 = feval(rbf4,epsilon,DM_eval);
%basisfunction5 = feval(rbf5,epsilon,DM_eval);

% Plot basis function
figure
hold on
plot(evalpoints, basisfunction1,'r-');
plot(evalpoints, basisfunction2,'g-.');
%plot(evalpoints, basisfunction3,'b--');
%plot(evalpoints, basisfunction4,'c-');
%plot(evalpoints, basisfunction5,'m:');
set(gca,'Fontsize',14)
xlabel('x','FontSize',14);
ylabel('y','FontSize',14,'Rotation',0);
%ylim([0 1])
%legend('C^0','C^2', 'C^4', 'C^6') % Wendland
%legend('C^0','C^2', 'C^4', 'C^6') % Wu
%legend('l=7/2','l=5', 'l=15/2', 'l=12') % Gneiting1
%legend('C^2', 'C^4', 'C^6') % Gneiting2
%legend('s=1','s=2', 's=3', 's=4', 's=5') % Euclid
caption = sprintf('RBF basis functions centered at origin.');
%title(caption)
hold off


