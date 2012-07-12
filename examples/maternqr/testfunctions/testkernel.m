  function out = testkernel(x)
  counter = 1;
  y = .5;
  out = zeros(length(x),1);
  for i = 1:length(x)
        if (0 < x(i) )&& (x(i) <= y)
            out(counter) = 
        elseif ((y <= x(i)) && (x(i) < 1))
            out(counter) = (-0.0696244 + exp(2.*x(i)).* (0.0909832 - 0.0407803.*x(i)) - 0.301328.*x(i));
        end
        counter = counter + 1;
  end
        