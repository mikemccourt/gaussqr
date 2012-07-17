function out = testfunction5(x)


out = zeros(length(x),1);
a = 4*pi;
b=1;
for i = 1:length(x)
    if ((x(i) >= 0) && (x(i)<= .5))
        out(i) = (.5/a) - (.5.*cos(a.*(x(i)).^2)/a);
    elseif ((x(i) <= 1) && (x(i)>.5))
        out(i) = ((0.5*cos(0.25.*a.*(1-2.*b).^2))/a)- ((0.5.*cos(a.*(b - 1.*x(i)).^2))/a) + (sin(0.125.*a).^2/a);
    end
end

end