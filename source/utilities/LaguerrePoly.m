% At some point, see if there's a way to not do this loop

function lp = LaguerrePoly(n, x, a)
lp = zeros(size(x));
avoid_zero = 1e-100;
for i=0:n
    lp = lp + (-1) ^ i * exp( ...
        gammaln(n + a + 1) + i * log(x + avoid_zero) - ...
        gammaln(n - i + 1) - gammaln(a + i + 1) - gammaln(i + 1) ...
    );
end
end