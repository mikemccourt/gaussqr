function out = testfunction4(x)


out = zeros(length(x),1);

for i = 1:length(x)
    if ((x(i) >= 0) && (x(i)<= .5))
        out(i) = sin(3.*x(i));
    elseif ((x(i) <= 1) && (x(i)>.5))
        out(i) = -sin(3.*(x(i)-1));
    end
end

end