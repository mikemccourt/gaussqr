% ex3
% This example shows the need for an asymptotic approximation of the
% Hermite polynomials as n->inf.  There may also be a need as x->inf, but
% we're not there yet.

figure

% For smallish n, polyval and asymptotic are pals
n=59;
xx=linspace(-2*sqrt(2*n),2*sqrt(2*n),100)';
hold on
plot(xx/sqrt(2*n),log(abs(HermiteAppx(n,xx))),'LineWidth',2)
plot(xx/sqrt(2*n),log(abs(polyval(HermitePoly(n),xx))),':r','LineWidth',2)
xlabel('x/sqrt(2n)')
ylabel('log(abs(H))')
title(sprintf('Hermite Polynomial of Order %d',n))
legend('Asymptotic','polyval','Location','North')
hold off

pause
clf reset

% For larger n, polyval will loss accuracy in a limited range
n=109;
xx=linspace(-2*sqrt(2*n),2*sqrt(2*n),100)';
hold on
plot(xx/sqrt(2*n),log(abs(HermiteAppx(n,xx))),'LineWidth',2)
plot(xx/sqrt(2*n),log(abs(polyval(HermitePoly(n),xx))),':r','LineWidth',2)
xlabel('x/sqrt(2n)')
ylabel('log(abs(H))')
title(sprintf('Hermite Polynomial of Order %d',n))
legend('Asymptotic','polyval','Location','North')
hold off

pause
clf reset

% For larger x, polyval will reach overflow
n=209;
xx=linspace(-2*sqrt(2*n),2*sqrt(2*n),100)';
hold on
plot(xx/sqrt(2*n),HermiteAppx(n,xx,1),'LineWidth',2)
plot(xx/sqrt(2*n),log(abs(polyval(HermitePoly(n),xx))),':r','LineWidth',2)
xlabel('x/sqrt(2n)')
ylabel('log(abs(H))')
title(sprintf('Hermite Polynomial of Order %d',n))
legend('Asymptotic','polyval','Location','North')
hold off

pause
clf reset

% For large enough n, all Hermite coefficients suffer overflow
n=309;
xx=linspace(-2*sqrt(2*n),2*sqrt(2*n),100)';
hold on
plot(xx/sqrt(2*n),HermiteAppx(n,xx,1),'LineWidth',2)
plot(xx/sqrt(2*n),log(abs(polyval(HermitePoly(n),xx))),':r','LineWidth',2)
xlabel('x/sqrt(2n)')
ylabel('log(abs(H))')
title(sprintf('Hermite Polynomial of Order %d',n))
legend('Asymptotic','polyval','Location','North')
hold off

Min_and_Max_coef_for_Hermite_309=[min(HermitePoly(309)),max(HermitePoly(309))]