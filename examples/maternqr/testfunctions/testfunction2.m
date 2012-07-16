function out = testfunction2(x)
%Comments!!!
%make lots
%of
%comments!

out = zeros(length(x),1);

for i = 1:length(x)
    if (0 <= x(i) )&& (x(i) < .1)
        out(i) = x(i).^2;
    elseif ((x(i)>= .1)&&(x(i)<= .9))
        out(i) = x(i).^3 -x(i).^2 + 3;
    elseif ((x(i)>.9)&&(x(i) <= 1))
        out(i) = (x(i) - 1).^2;
    else
        x(i) = 0;
    end
end

end