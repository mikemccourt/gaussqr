%errc = erroref(Phi,M)
%Evaluates the error for the eigenfunction
%
%Inputs : Phi - The object get form HSeigsolve
%         M   - Methods that will be used to compute error.
%
%Outputs %errc - The error for eigenfunction

function errc = erroref(Phi,M)
% compute eigenfunctions' realtive error
 C = Phi.coefs;
 sc = size(C);
 j = sc(2); %how many columns
 rl = sc(1); %how many rows   
 errc = zeros(rl,1);
   switch M
       case 1 %abs(normone(realfunction)-normone(apporximatefunction))
           for i =1: j
               errc(i) = abs( 2*sqrt(2)/pi - quad(@(x) Phi.eigfuncEval(x',i)',0,1));
           end
       case 2 %normone(realfunction - apporximatefunction)
           for i =1: j
               errc(i) = quad(@(x) abs(sqrt(2)*sin(pi*i*x') - Phi.eigfuncEval(x',i)),0,1);
           end
       case 3 %boundary value
            for i =1: j
               errc(i) = abs(Phi.eigfuncEval(0,i))/sqrt(2); %max value of function ff should be sqrt(2)
             end
        case 4 %collection points
             for i =1: j
                errc(i) =sum(abs(sqrt(2)*sin(pi*i*(0:0.002:1)')-Phi.eigfuncEval((0:.002:1)',i)))/500*pi/(2*sqrt(2)); 
             end
        otherwise error('Unacceptable method=%e',M);
    end