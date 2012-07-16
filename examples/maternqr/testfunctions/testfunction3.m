function out = testfunction3(x)


out = zeros(length(x),1);

for i = 1:length(x)
    if (x(i) ~= (0))
        out(i) = (x(i).^2).*(sin(1/x(i)));
    elseif (x(i) == 0)
        out(i) = 0;
    end
end

end