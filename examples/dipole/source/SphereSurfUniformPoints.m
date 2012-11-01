function [vertici, triangoli] = SphereSurfUniformPoints( r, maxlivello )
% Genera una superficie composta da triangoli che approssima una sfera di
% raggio r tramite suddivisioni successive dei lati di un un icosaedro.
% maxlivello fissa il numero di divisioni successive (meglio non andare 
% oltre 4).

% Coordinate dei dodici vertici dell'icosaedro (riferimento a sfera
% unitaria)
tau = 0.8506508084; % t=(1+sqrt(5))/2, tau=t/sqrt(1+t^2)
one = 0.5257311121; % one=1/sqrt(1+t^2)

vertici( 1,:) = [  tau,  one,    0 ]; % ZA
vertici( 2,:) = [ -tau,  one,    0 ]; % ZB
vertici( 3,:) = [ -tau, -one,    0 ]; % ZC
vertici( 4,:) = [  tau, -one,    0 ]; % ZD
vertici( 5,:) = [  one,   0 ,  tau ]; % YA
vertici( 6,:) = [  one,   0 , -tau ]; % YB
vertici( 7,:) = [ -one,   0 , -tau ]; % YC
vertici( 8,:) = [ -one,   0 ,  tau ]; % YD
vertici( 9,:) = [   0 ,  tau,  one ]; % XA
vertici(10,:) = [   0 , -tau,  one ]; % XB
vertici(11,:) = [   0 , -tau, -one ]; % XC
vertici(12,:) = [   0 ,  tau, -one ]; % XD

% Struttura dell'iscosaedro unitario
triangoli = [  5,  8,  9 ;
               5, 10,  8 ;
               6, 12,  7 ;
               6,  7, 11 ;
               1,  4,  5 ;
               1,  6,  4 ;
               3,  2,  8 ;
               3,  7,  2 ;
               9, 12,  1 ;
               9,  2, 12 ;
               10,  4, 11 ;
               10, 11,  3 ;
               9,  1,  5 ;
               12,  6,  1 ;
               5,  4, 10 ;
               6, 11,  4 ;
               8,  2,  9 ;
               7, 12,  2 ;
               8, 10,  3 ;
               7,  3, 11 ];

% Affinamento della griglia tramite suddivisioni
for level = 1:maxlivello,
	% Suddivisione di ogni triangolo
    [vertici, triangoli] = affinamento(vertici, triangoli);  
    % Proiezione dei punti sulla superficie di una sfera di raggio r
	vertici = proiez_sfera(vertici,r); 
end
triangoli = triangoli(:,[1 3 2]); % Ordine dei vertici antiorario guardando 
                                  % dall'esterno della superficie
end

function [vertici, triangoli] = affinamento(vertici, triangoli)
% Infittisce la griglia introducendo nuovi punti in corrispondenza dei
% punti medi dei lati dei triangoli già presenti:
%
%        B
%        /\
%      a/__\b
%      /\  /\
%     /__\/__\
%    A   c    C
% Nuovi punti medi:
% a = (A+B)/2
% b = (B+C)/2
% c = (C+A)/2
%
% Struttura triangoli:
% [A,a,c]
% [a,B,b]
% [c,b,C]
% [a,b,c]

% Inizializzazione nuove matrici dei vertici e dei triangoli
Nface = size(triangoli,1);
F2 = zeros(Nface*4,3);

for f = 1:Nface,
	% Lettura indici dei vertici del triangolo
    NA = triangoli(f,1);
    NB = triangoli(f,2);
    NC = triangoli(f,3);   
    
    % Lettura coordinate dei vertici del triangolo
    A = vertici(NA,:);
    B = vertici(NB,:);
    C = vertici(NC,:);
    
    % Nuovi vertici (punti medi)
    a = (A + B) ./ 2;
    b = (B + C) ./ 2;
    c = (C + A) ./ 2;
    
    % Immagazzinamento nuovi vertici e controllo esistenza
    [vertici, Na] = trova_vertice(vertici,a);
    [vertici, Nb] = trova_vertice(vertici,b);
    [vertici, Nc] = trova_vertice(vertici,c); 
    
    % Creazione nuovi triangoli
    F2(f*4-3,:) = [ NA, Na, Nc ];
    F2(f*4-2,:) = [ Na, NB, Nb ];
    F2(f*4-1,:) = [ Nc, Nb, NC ];
    F2(f*4-0,:) = [ Na, Nb, Nc ];
end

% Nuova matrice dei triangoli
triangoli = F2;
end

function [vertici, N] = trova_vertice( vertici, vertice )
    Vn = size(vertici,1);
    Va = repmat(vertice,Vn,1);
    Vesist = find( vertici(:,1) == Va(:,1) & ...
                   vertici(:,2) == Va(:,2) & ...
                   vertici(:,3) == Va(:,3) );
    if Vesist,
        if size(Vesist) == [1,1],
            N = Vesist;
        else
            error('Vertici replicati');
        end
    else
        vertici(end+1,:) = vertice;
        N = size(vertici,1);
    end
return
end

function [vertici] = proiez_sfera( v, r, c )
% Proietta i punti X,Y,Z su una sfera di raggio r
%
% Utilizzo:
%  vertici = proiez_sfera( v, r, c )
%
% in cui:
% v è la matrice dei vertici, N x 3 (XYZ)
% r è il raggio della sfera, 1 x 1 (di default 1)
% c è il centro della sfera, 1 x 3 (di default 0,0,0)
%
% XYZ vengono convertiti in coordinate sferiche e i rispettivi raggi
% aggiunstati in accordo con il valore r indicato dal centro c verso XYZ
% (attraverso le coordinate phi & theta).
%
% I vertici vengono forniti in coordinate cartesiane

if ~exist('v','var'),
    error('proiez_sfera: Nessun vertice in ingresso (X,Y,Z)\n');
end

X = v(:,1);
Y = v(:,2);
Z = v(:,3);

if ~exist('c','var'), % se il centro non è indicato, allora viene posto 
                      % nell'origine
    xo = 0;
    yo = 0;
    zo = 0;
else
    xo = c(1);
    yo = c(2);
    zo = c(3);
end

if ~exist('r','var'), r = 1; end % se il raggio non è indicato, allora 
                                 % viene posto r=1
 
% Converte le coordinate cartesiane X,Y,Z in sferiche theta e phi
theta = atan2( (Y-yo), (X-xo) );
phi   = atan2( sqrt( (X-xo).^2 + (Y-yo).^2 ), (Z-zo) );

% Ricalcola X,Y,Z (proiezioni dei vertici dell'icosaedro sulla sfera)
R = ones(size(phi)) * r;
x = R .* sin(phi) .* cos(theta);
y = R .* sin(phi) .* sin(theta);
z = R .* cos(phi);

% Nuova matrice dei vertici
vertici = [x y z];
end

