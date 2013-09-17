%% Validation/regression tests
%
% This is NOT a variety of regression/least squares tool. In fact, it
% has nothing to do with those capabilities of MATLAB. Instead, 
% this document provides tests for each method in the HPF world, so
% that any bugs in future changes/enhancements to HPF will be efficiently
% caught. Each method in HPF is included, with a set of tests chosen
% specifically for that method to exercise the numerical precision and
% behavior of the method. Special cases are tested, such as the behavior of
% each function on empty arrays, inf, nan, etc. to conform to the expected
% behavior in MATLAB.
%
% As a user, you can publish this file yourself on your own system, to
% then verify that all tests return the expected values.
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com

format long g
%%
% Save the current state to restore it later on
OriginalDecimalBase = DefaultDecimalBase;
OriginalNDig = DefaultNumberOfDigits;


% Set the default digits to [64 4]
DefaultNumberOfDigits 64 4

%% abs: absolute value
% empty begets empty
abs(hpf([]))
%%
% Goals: 12.34567890 0 1.234 inf nan
X = [hpf('-12.34567890'), hpf(0), hpf('1.234'), hpf(-inf), hpf(nan)];
abs(X)

%% acos: inverse cosine
% empty begets empty
acos(hpf([]))
%%
% Goals: 0, 1/4, 1/2, 2/3, 1, nan nan nan
X = [hpf(1), sqrt(hpf(2))/2, hpf(0), hpf(-0.5), hpf(-1), hpf(-inf), hpf(inf), hpf(nan)];
acos(X)/hpf('pi')

%%
% A set of random tests.
% The goal is a max log10 difference of roughly -15
x = hpf(2*rand(1,25) - 1);
log10(max(abs(acos(x) - acos(double(x)))))

%%
% Sym comparison:
%
% vpa(acos(sym('1/4')),100)
%
% ans = 1.318116071652817965745664254646040469846390966590714716853548517413333142662083276902268670443043932
acos(hpf('0.25',100))

%% acosd:
% empty begets empty
acosd(hpf([]))
%%
% Goals: 0, 1/4, 1/2, 2/3, 1, nan nan nan
X = [hpf(1), sqrt(hpf(2))/2, hpf(0), hpf(-0.5), hpf(-1), hpf(-inf), hpf(inf), hpf(nan)];
acosd(X)/180

%%
% A set of random tests.
% The goal is a max log10 difference of roughly -13
x = hpf(2*rand(1,10) - 1,100);
log10(max(abs(acosd(x) - acosd(double(x)))))

%% acosh:
% empty begets empty
acosh(hpf([]))
%%
% Goals: [0, inf, NaN]
acosh(hpf([1 inf NaN],50))

%%
% A set of random tests.
% The goal is a max log10 difference of roughly -15
x = hpf(exp(2*rand(1,10)));
log10(max(abs(double(acosh(x)) - acosh(double(x)))))

%%
% Sym comparison:
%
% vpa(acosh(sym(17)),100)
%
% ans = 3.525494348078172100930437299919169236112641313046541643013182434613508736888104351334827567641024171
acosh(hpf('17',100))


%% acot:
% empty begets empty
acot(hpf([]))

%%
% Goals: [pi/2, 0, 0, NaN]
acot(hpf([0 -inf inf NaN],50))

%%
% A set of random tests.
% The goal is a max log10 difference of roughly -15
x = hpf(10*rand(1,10)-5,50);
log10(max(abs(double(acot(x)) - acot(double(x)))))

%%
% Sym comparison:
%
% vpa(acot(sym(17)),100)
%
% ans = 0.05875582271572269290218330245177178547665041810044464614350256296744422042319255800245502374523156074
acot(hpf('17',100))

%% acotd:




%% acsc:
% empty begets empty
acsc(hpf([]))

%%
% Goals: [-pi/2, pi/2, 0, 0, NaN]
acsc(hpf([-1 1 -inf inf NaN],50))

%%
% A set of random tests.
% The goal is a max log10 difference of roughly -15
x = hpf(5*rand(1,20).*(2*round(rand(1,20))-1),30);
log10(max(abs(double(acsc(x)) - acsc(double(x)))))

%%
% Sym comparison:
%
% vpa(acsc(sym(17)),100)
%
% ans = 0.05885750594708123021242193548567836102172841359730737506627786418916451737407051526190724208378293111
acsc(hpf('17',100))

%% acscd:




%% adjustdecimalbase:




%% asec
% empty begets empty
asec(hpf([]))

%%
% Goals: [pi, 0, pi/2, pi/2, NaN]
asec(hpf([-1 1 -inf inf NaN],50))

%%
% A set of random tests.
% The goal is a max log10 difference of roughly -15
x = hpf(10*rand(1,10)-5,30);
log10(max(abs(double(asec(x)) - asec(double(x)))))

%%
% Sym comparison:
%
% vpa(asec(sym(17)),100)
%
% ans = 1.511938820847815389018899756154073081076856286090245535421194431964743685769033984052110170587275603
asec(hpf('17',100))

%% asecd:







%% asin:
% Goals: [], NaN, NaN, NaN ,NaN
asind(hpf([]))
asind(hpf([2 -inf inf NaN]))

%%
% Goals: [0 1/6 -1/2]
asin(hpf([0 0.5 -1],[200,3]))/hpf('pi',[200 3])

%%
% A set of random tests.
% The goal is a max log10 difference of roughly -15
x = hpf(2*rand(1,25) - 1,[50,2]);
log10(max(abs(asin(x) - asin(double(x)))))

%%
% Sym comparison:
%
% vpa(sin(sym('3')),100)
%
% ans = 0.1411200080598672221007448028081102798469332642522655841518826412324220099670144719112821728534498638
sin(hpf('3',100))

%% asind:
% Goals: [], NaN, NaN ,NaN
asind(hpf([]))
asind(hpf([-inf inf NaN]))

%%
% Goals: [0 30 -90]
asind(hpf([0 0.5 -1],[200,3]))

%%
% A set of random tests.
% The goal is a max log10 difference of roughly -13
x = hpf(2*rand(1,25) - 1,[30 2]);
log10(max(abs(asind(x) - asind(double(x)))))

%% asinh:
% empty begets empty
asinh(hpf([]))
%%
% Goals: [-inf, 0, inf, NaN]
asinh(hpf([-inf 0 inf NaN],50))

%%
% A set of random tests.
% The goal is a max log10 difference of roughly -15
x = hpf(20*rand(1,10) - 10);
log10(max(abs(double(asinh(x)) - asinh(double(x)))))

%%
% Sym comparison:
%
% vpa(asinh(sym(17)),100)
%
% ans = 3.527224456199965806371844321147769737181075325582383725627847771483069723308938434014739273560152791
asinh(hpf('17',100))

%% atan:
% empty begets empty
atan(hpf([]))
%%
% Goals: [-1, 1, NaN]
atan(hpf([-inf inf NaN],50))./(hpf('pi',50)./2)

%%
% Goals: [0 -0.25]
atan(hpf([0 -1],[200,3]))/hpf('pi',[200 3])

%%
% Goal: 1
atan(2 - sqrt(hpf(3,100)))*12/hpf('pi',100)

%%
% A set of random tests.
% The goal is a max log10 difference of roughly -15
x = hpf(20*rand(1,10) - 10);
log10(max(abs(double(atan(x)) - atan(double(x)))))

%%
% Sym comparison:
%
% vpa(atan(sym(17)),100)
%
% ans = 1.512040504079173926329138389187979656621934281587108264343969733186463982719911941311562388925826973
atan(hpf('17',100))

%% atand:





%%
% A set of random tests.
% The goal is a max log10 difference of roughly -13
x = hpf(20*rand(1,25) - 10);
log10(max(abs(atand(x) - atand(double(x)))))

%% atanh:
% empty begets empty
atanh(hpf([]))
%%
% Goals: [-inf, 0, inf, NaN]
atanh(hpf([-1 0 1 NaN],50))

%%
% A set of random tests.
% The goal is a max log10 difference of roughly -15
x = hpf(2*rand(1,10) - 1);
log10(max(abs(double(atanh(x)) - atanh(double(x)))))

%%
% Sym comparison:
%
% vpa(atanh(sym(0.5)),100)
%
% ans = 0.5493061443340548456976226184612628523237452789113747258673471668187471466093044834368078774068660444
atanh(hpf('0.5',100))


%% augmentdigits:




%% ceil:




%% cos: cosine function
% empty begets empty
cos(hpf([]))

%%
% Goals [1 -1 sqrt(2)/2 1/2 nan nan nan]
X= [hpf(0), 3*hpf('pi',200), hpf('pi')/4, hpf('pi')/3, hpf(-inf) hpf(inf) hpf(nan)];
cos(X)
cos(double(X))

%%
% A set of random tests.
% The goal is a max log10 difference of roughly -15
x = hpf(2*pi*rand(1,25));
log10(max(abs(cos(x) - cos(double(x)))))

%%
% Sym comparison:
%
% vpa(cos(sym('3')),100)
%
% ans = -0.9899924966004454572715727947312613023936790966155883288140859329283291975131332204282944793556926022
cos(hpf('3',100))

%% cosd:








%% cosh:
% empty begets empty
sinh(hpf([]))

%%
% Goals: [1 inf inf NaN
cosh(hpf([0 -inf inf NaN]))

%%
% Fixed tests: Goal is -15.288164505228
x = [-3 -2 -1 1 2 3];
log10(max(abs(double(cosh(hpf(x,[500 4])) - cosh(x)))))

%%
% Sym comparison:
%
% vpa(cosh(sym(1)),100)
%
% ans = 1.543080634815243778477905620757061682601529112365863704737402214710769063049223698964264726435543036
cosh(hpf(1,100))

%% coth:
% empty begets empty
sinh(hpf([]))

%%
% Goals: [Inf -1 1 NaN
coth(hpf([0 -inf inf NaN]))

%%
% Fixed tests: Goal is -15.8016545031234
x = [-3 -2 -1 1 2 3];
log10(max(abs(double(coth(hpf(x,[500 4])) - coth(x)))))

%%
% Sym comparison:
%
% vpa(coth(sym(1)),100)
%
% ans = 1.313035285499331303636161246930847832912013941240452655543152967567084270461874382674679241480856303
coth(hpf(1,100))

%% csc





%% cscd





%% csch:




%% cubrt
% Goals: [0, 0, -1.1, 8.37, inf, -inf, NaN]
cubrt(exp(hpf('3',500))) - hpf('e',500)
cubrt(hpf('-1.331',1000))
cubrt(hpf('586.376253',250))
cubrt(hpf([inf -inf NaN]))

%% cumprod




%% cumsum





%% DefaultNumberOfDigits: sets the default working precision of HPF






%% DefaultDecimalBase: sets the default working base for HPF
% Goal
sqrt(pi*exp(1) + 2) - 4/3

%%
% Decimal bases of 1:6, should all yield the same result
DefaultNumberOfDigits 50 1
DefaultDecimalBase 1
sqrt(hpf('pi')*hpf('e') + 2) - hpf(4)/3

DefaultDecimalBase 2
sqrt(hpf('pi')*hpf('e') + 2) - hpf(4)/3

DefaultDecimalBase 3
sqrt(hpf('pi')*hpf('e') + 2) - hpf(4)/3

DefaultDecimalBase 4
sqrt(hpf('pi')*hpf('e') + 2) - hpf(4)/3

DefaultDecimalBase 5
sqrt(hpf('pi')*hpf('e') + 2) - hpf(4)/3

DefaultDecimalBase 6
sqrt(hpf('pi')*hpf('e') + 2) - hpf(4)/3

%%
% A time test, with 35000 digits.
% Time for a multiply in HPF should be roughly inversely quadratic,
% but adds and subtracts are not so strongly a function of the decimal
% base, although at worst linear in the number of digits, and divides
% (and other ops like sqrt) will be a mixture of adds and multiplies.

%%
% Time tests for a pure add
DefaultNumberOfDigits 35000 0
DefaultDecimalBase 1
timeit(@() hpf('pi') + hpf('e'))

DefaultDecimalBase 2
timeit(@() hpf('pi') + hpf('e'))

DefaultDecimalBase 3
timeit(@() hpf('pi') + hpf('e'))

DefaultDecimalBase 4
timeit(@() hpf('pi') + hpf('e'))

DefaultDecimalBase 5
timeit(@() hpf('pi') + hpf('e'))

DefaultDecimalBase 6
timeit(@() hpf('pi') + hpf('e'))

%%
% Time tests for a pure multiply
DefaultDecimalBase 1
timeit(@() hpf('pi').*hpf('e'))

DefaultDecimalBase 2
timeit(@() hpf('pi').*hpf('e'))

DefaultDecimalBase 3
timeit(@() hpf('pi').*hpf('e'))

DefaultDecimalBase 4
timeit(@() hpf('pi').*hpf('e'))

DefaultDecimalBase 5
timeit(@() hpf('pi').*hpf('e'))

DefaultDecimalBase 6
timeit(@() hpf('pi').*hpf('e'))

%%
% Time tests for a pure divide
DefaultDecimalBase 1
tic, hpf('pi')./hpf('e');toc

DefaultDecimalBase 2
tic, hpf('pi')./hpf('e');toc

DefaultDecimalBase 3
tic, hpf('pi')./hpf('e');toc

DefaultDecimalBase 4
tic, hpf('pi')./hpf('e');toc

DefaultDecimalBase 5
tic, hpf('pi')./hpf('e');toc

DefaultDecimalBase 6
tic, hpf('pi')./hpf('e');toc

%%
% Time tests for a mixed operation
DefaultDecimalBase 1
tic,sqrt(hpf('pi')*hpf('e') + 2) - hpf(4)/3;toc

DefaultDecimalBase 2
tic,sqrt(hpf('pi')*hpf('e') + 2) - hpf(4)/3;toc

DefaultDecimalBase 3
tic,sqrt(hpf('pi')*hpf('e') + 2) - hpf(4)/3;toc

DefaultDecimalBase 4
tic,sqrt(hpf('pi')*hpf('e') + 2) - hpf(4)/3;toc

DefaultDecimalBase 5
tic,sqrt(hpf('pi')*hpf('e') + 2) - hpf(4)/3;toc

DefaultDecimalBase 6
tic,sqrt(hpf('pi')*hpf('e') + 2) - hpf(4)/3;toc

%% disp:




%% display:




%% double:




%% eps:






%% eq:
DefaultNumberOfDigits 64 4
X = hpf([-inf inf NaN 3   2.5 -1.37 5 6.2999999999999 6.4]);
Y = hpf([inf inf  inf NaN -1   0    2 6.3999999999999 6.4]);
X == Y
%%
% Expected result
double(X) == double(Y)


%% erf
% empty begets empty
erf(hpf([]))

%%
% Goals: [-1 0 1 NaN]
erf(hpf([-inf 0 inf NaN]))

%%
% Fixed tests: Goal is -16.2705600471626
x = -5:2:5;
log10(max(abs(double(erf(hpf(x,[20 4])) - erf(x)))))

%%
% Sym comparison:
%
% vpa(erf(sym(9)),100)
%
% ans = 0.9999999999999999999999999999999999995862968253486189761946096532637475404289808014052037540194090098
erf(hpf('9',100))

%% exp:
% empty begets empty
exp(hpf([]))

%%
% Goals: [0,1,exp(1),inf,nan,10]
X = [hpf([-inf 0 1 inf nan]),hpf('ln10')];
exp(X)

%%
% Fixed tests: Goal is -14.4576292646912
x = -5:2:5;
log10(max(abs(double(exp(hpf(x,[20 4])) - exp(x)))))

%%
% Sym comparison:
%
% vpa(exp(sym('5/2')),100)
%
% ans = 12.18249396070347343807017595116796618318276779006316131156039834183818512614331441006025552300629579
exp(hpf('2.5',100))

%% factorial:
% factorial(sym(75))
%
% ans = 24809140811395398091946477116594033660926243886570122837795894512655842677572867409443815424000000000000000000
factorial(hpf(75,200))

%% fix:
% empty begets empty
fix(hpf([]))

%%
% Goals: [-inf inf NaN]
fix(hpf([-inf inf NaN],23))

%%
% Goal
Xd = (-10:10)*pi/10;
fix(Xd)
%%
fix(hpf(Xd,[20,0]))
%%
% Should be entirely true:
double(fix(hpf(Xd,[20,0]))) == fix(Xd)

%% floor: round towards -inf
% Empty begets empty
floor(hpf([]))

%%
% Expected results: [3,-314160 12345678901 -1234567890123456789 0 -inf NaN]
X = [hpf('pi'), -1e5*hpf('pi'), hpf('12345678901'), ...
  hpf('-123456789012345678901'), hpf([0, -inf, nan])]
floor(X)

%% fractionalpart: 
% Empty begets empty
fractionalpart(hpf([]))

%%
% Goals: [0.141590, -0.345 0 0 NaN NaN NaN]
DefaultDecimalBase 6
fractionalpart(hpf('pi',[6 0]))
fractionalpart(hpf('-12.345'))
fractionalpart(hpf([0 10000 -inf inf NaN]))

%% ge: 
X = hpf([-inf inf NaN 3   2.5 -1.37 5 6.2999999999999 6.4]);
Y = hpf([inf inf  inf NaN -1   0    2 6.3999999999999 6.4]);
X >= Y
%%
% Expected result
double(X) >= double(Y);

%%
% Random tests, Goal: 1
x = hpf(1 + (rand(1,10) - 0.5)/100000);
y = hpf(1 + (rand(1,10) - 0.5)/100000);
all((x >= y) == (double(x) >= double(y)))


%% gt: 
X = hpf([-inf inf NaN 3   2.5 -1.37 5 6.2999999999999 6.4]);
Y = hpf([inf inf  inf NaN -1   0    2 6.3999999999999 6.4]);
X > Y
%%
% Expected result
double(X) > double(Y)

%%
% Random tests, Goal: 1
x = hpf(1 + (rand(1,10) - 0.5)/100000);
y = hpf(1 + (rand(1,10) - 0.5)/100000);
all((x > y) == (double(x) > double(y)))

%% hpf: 
% empty begets empty
hpf([])

%%
% Goal: 0
hpf

%%
% Goals: [0 1 -inf inf NaN]
hpf([0 1 -inf inf NaN])

%%
% Goals: [1.23, 1.229999999999999982236431605997495353221893310546875000]
hpf('1.23')
hpf(1.23)

%%
% Goals: [123456.7, 12345.67, 12.34567, 1.234567e10000, -12.1212121]
hpf('1.234567e5')
hpf('0.1234567e5')
hpf('1234.567e-2')
hpf('1234567e10000')
hpf('-12.1212121e0')

%%
% Goal: 1.2345678901234567890123e29
hpf('123456789012345678901234567890',23)

%% int16




%% int32




%% int64





%% int8







%% isfinite: 
% empty begets empty
isfinite(hpf([]))

%%
% goals: [1 1 1 0 0 0]
isfinite(hpf([0 pi -4 NaN inf -inf]))

%% isinf: 
% empty begets empty
isinf(hpf([]))

%%
% goals: [0 0 0 0 1 1]
isinf(hpf([0 pi -4 NaN inf -inf]))

%% isnan:
% empty begets empty
isnan(hpf([]))

%%
% goals: [0 0 0 1 0 0]
isnan(hpf([0 pi -4 NaN inf -inf]))

%% isnumeric:
% Goal: 1 (true)
isnumeric(hpf)

%% le: <=, less than or equal to
X = hpf([-inf inf NaN 3   2.5 -1.37 5 6.2999999999999 6.4]);
Y = hpf([inf inf  inf NaN -1   0    2 6.3999999999999 6.4]);
X <= Y
%%
% Expected result
double(X) <= double(Y)







%% linspace




%% log: 







%% log10: 
% Goals: [], -inf, inf, NaN
log10(hpf([]))
log10(hpf(0))
log10(hpf(inf))
log10(hpf(NaN))

%%
% Goals: [50, 100000000000000000000000000000000000000000000000001]
DefaultNumberOfDigits 200 5
X = hpf('100000000000000000000000000000000000000000000000000');
log10(X)
10.^(log10(X+1))

%%
% Sym comparison:
%
% vpa(log10(sym('17')),100)
%
% ans = 1.23044892137827392854016989432833703000756737842504639738036848234469406225711818579568467009846514
log10(hpf('17',100))

%% log2: 
% Goals: [], -inf, inf, NaN
log2(hpf([]))
log2(hpf(0))
log2(hpf(inf))
log2(hpf(NaN))

%%
% Goals: [19520]
X = hpf('1048576',[50 5]).^976;
log2(X)

%%
% Sym comparison:
%
% vpa(log2(sym('17')),100)
%
% ans = 4.087462841250339408254066010810404354011267282344820688126609064386696509047382068297343151843684273
log2(hpf('17',100))

%% lt: <, less than
X = hpf([-inf inf NaN 3   2.5 -1.37 5 6.2999999999999 6.4]);
Y = hpf([inf inf  inf NaN -1   0    2 6.3999999999999 6.4]);
X < Y
%%
% Expected result
double(X) < double(Y)



%% mantissa:




%% max




%% min




%% minus:





%% mod: 





%% mpower: 






%% mrdivide: 





%% mtimes: 






%% ne: 
X = hpf([-inf inf NaN 3   2.5 -1.37 5 6.2999999999999 6.4]);
Y = hpf([inf inf  inf NaN -1   0    2 6.3999999999999 6.4]);
X ~= Y
%%
% Expected result
double(X) ~= double(Y)

%% plus: 
% empty begets empty, as an empty HPF result
[] + hpf(pi)

%%
% Goals: [NaN, inf, NaN, NaN, 1.5 0], 
X = hpf([-inf inf NaN 3   2.5 pi]);
Y = hpf([inf inf  inf NaN -1  -pi]);
X + Y
%%
% Goal: 1.234567, with NumberOfDigts == [6 0]
hpf('1.11111111111111111111',[20 0]) + hpf('0.123456',[6 0]) 
ans.NumberOfDigits

%% nthroot



%% plus



%% power: 




%% prod




%% rdivide: 






%% reciprocal: Computes the scalar inverse of a number, 1./F
% empty begets empty
reciprocal(hpf([]))

%%
% Goals: [inf, NaN, 0, 0]
reciprocal(hpf([0, NaN, -inf, inf]))

%%
% a set of tests on reciprocals of some integers. Goal is -16.255619765855
log10(max(abs(double(reciprocal(hpf(1:50,20))) - 1./(1:50))))

%%
% A set of reciprocals of some random numbers. Goal is roughly -14
x = sort(rand(1,20) - 1);
log10(max(abs(double(reciprocal(hpf(x,[30 2]))) - 1./x)))

%% round: Round to the nearest integer



%% roundn: Round to the nearest indicated power of 10



%% sec
% empty begets empty
sec(hpf([]))

%%
% Goals: [inf inf 1 NaN NaN NaN]
sec(-hpf('pi',40)/2)
sec(hpf('pi',40)/2)
sec(hpf([0 -inf inf NaN]))

%%
% Fixed tests: Goal is  -15.083071513751
x = 1:10;
log10(max(abs(double(sec(hpf(x,[20 4])) - sec(x)))))

%%
% Sym comparison:
%
% vpa(sec(sym(17)),100)
%
% ans = -3.634205076449851045700145553853416971812501547277707974413306626204205786111854275175875759446174066
sec(hpf('17',100))

%% secd



%% sech: 
% empty begets empty
sech(hpf([]))

%%
% Goals: [0 1 0 NaN]
sech(hpf([-inf 0 inf NaN]))

%%
% Fixed tests: Goal is  -15.083071513751
x = 1:10;
log10(max(abs(double(sec(hpf(x,[20 4])) - sec(x)))))

%%
% Sym comparison:
%
% vpa(sech(sym(17)),100)
%
% ans = 0.00000008279875437570319128353730868526605265319669601947623809070883019333038402979377756328623635837203043
sech(hpf('17',100))

%% sign: 




%% sin: 
% empty begets empty
sin(hpf([]))

%%
% Goals: [0 NaN NaN NaN]
sin(hpf([0 -inf inf NaN]))

%%
% Fixed tests: Goal is -16.3104935951442
x = 1:10;
log10(max(abs(double(sin(hpf(x,[20 4])) - sin(x)))))

%%
% special values, goals: 0 (to within 250 digits)
x = hpf('pi',250)./[6 4 3 2 1]
sin(x) - [0.5 sqrt(hpf(2,250))/2 sqrt(hpf(3,250))/2 1 0]

%%
% Sym comparison:
%
% vpa(sin(sym(17)),100)
%
% ans = -0.961397491879556857261636944869156098492067254058935985601562874569876633233567787053089546435476294
sin(hpf('17',100))

%% sind: 





%% single: 







%% sinh: 
% empty begets empty
sinh(hpf([]))

%%
% Goals: [0 -inf inf NaN
sinh(hpf([0 -inf inf NaN]))

%%
% Random tests: Goal is -15.1562753989016
x = [-3 -2 -1 1 2 3];
log10(max(abs(double(sinh(hpf(x,[500 4])) - sinh(x)))))

%%
% Test of a large value. Goal: -2.57653587596...e-869
x = sinh(hpf(1000,[1000 5]));
log(x) + log(hpf('2',[1000 5])) - 1000

%%
% Sym comparison:
%
% vpa(sinh(sym(1)),100)
%
% ans = 1.175201193643801456882381850595600815155717981334095870229565413013307567304323895607117452089623392
sinh(hpf(1,100))



%% sort: 





%% sqrt: 
% empty begets empty
sqrt(hpf([]))

%%
% Goals: [0 NaN inf NaN]
%
% A warning will be issued for sqrt(-1)
sqrt(hpf([0 -1 inf NaN]))

%%
% Fixed tests: Goal is -15.6638172422041
x = 2:10;
log10(max(abs(double(sqrt(hpf(x,[20 4])) - sqrt(x)))))

%%
% Sym comparison:
%
% vpa(sqrt(sym(23)),100)
%
% ans = 4.795831523312719541597438064162693919996707041904129346485309114448257235907464082492191446436918861
sqrt(hpf('23',100))

%% sum







%% tan: 
% empty begets empty
tan(hpf([]))

%%
% Goals: [0 NaN NaN NaN]
tan(hpf([0 -inf inf NaN]))

%%
% Goals: [-inf inf]
tan(hpf('pi',100)*[-.5 .5])

%%
% Fixed tests: Goal is -15.5108867853615
x = 1:10;
log10(max(abs(double(tan(hpf(x,[20 4])) - tan(x)))))

%%
% special values, goals: 0 (to within 1e-150
x = hpf('pi',150)./[6 4 3]
tan(x) - [sqrt(hpf(3,150))/3 1 sqrt(hpf(3,250))]

%%
% Sym comparison:
%
% vpa(tan(sym(17)),100)
%
% ans = 3.493915645474839978346364290913980748486600343707345139970378981729932323805208911327226624805744553
tan(hpf('17',100))

%% tand: 





%% tanh: 
% empty begets empty
tanh(hpf([]))

%%
% Goals: [0 -1 1 NaN]
tanh(hpf([0 -inf inf NaN]))

%%
% Fixed tests: Goal is -15.8487333751523
x = 0.1:.1:1;
log10(max(abs(double(tanh(hpf(x,[20 4])) - tanh(x)))))

%%
% Sym comparison:
%
% vpa(tanh(sym(17)),100)
%
% ans = 0.9999999999999965721831369159799423587827364437255846875211707191426378969848776206130640331425577278
tanh(hpf('17',100))

%% times: 



%% uint16



%% uint32




%% uint64




%% uint8





%% uminus: unary minus, -F




%% uplus: 




%% vpi: 




%%
% Restore the System state to whatever it was before this test was run.
DefaultDecimalBase(OriginalDecimalBase)
DefaultNumberOfDigits(OriginalNDig)

