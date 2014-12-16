% FrontCover
% This is an attempt to create an image for the front cover
% The painting this is based on is available on the GaussQR website
gqr_downloaddata('coverpainting.png')

% Read the image into the computer
X = imread('coverpainting.png','png');
[ny,nx,nc] = size(X);
N = nx*ny;

% Convert the image to data for manipulation
% We rescale the data to [0,1], which is not necessary, but whatever
Xs = double(X);

% Choose a kernel to work with
% It will have to be compactly supported for a problem this big
% This should probably be applied with spfun
rbfW0 = @(e,r) max(1-e*r,0).^2;
rbfW2 = @(e,r) max(1-e*r,0).^4.*(4*e*r+1);
rbfW4 = @(e,r) max(1-e*r,0).^6.*(35*(e*r).^2+18*e*r+3);
rbfW6 = @(e,r) max(1-e*r,0).^8.*(32*(e*r).^3 + 25*(e*r).^2+8*e*r+1);

% Choose kernel parameters
% ep is the shape parameter, mu is for ridge regression
ep = 100.01;
mu = 1e-8;

% Create data locations, which will just be evenly spaced points
x = pick2Dpoints([0 0],[1 1],[nx;ny]);

% Form the necessary distance matrix now for use later
% This is complicated because this requires the use of a sparse matrix for
% a problem this big
% I am going to execute this computation in pieces so that I can monitor
% its progress
progress_tick = 1/2;
progress = 0;
h_waitbar = waitbar(progress,'Initializing distance matrix');

kstep = ceil(N*progress_tick);
krange = 1:kstep;
[Cidx,Cdist] = rangesearch(x,x(krange,:),1/ep);
nz_new = sum(cellfun(@length,Cidx));
% Allocated expected space for the i, j, s vectors which will be used to
% define the sparse matrix
% Then fill up those vectors with the 
ivec = zeros(1,nz_new/progress_tick);
jvec = zeros(1,nz_new/progress_tick);
svec = zeros(1,nz_new/progress_tick);
ivec(1:nz_new) = cell2mat(cellfun(@(x,k)k*ones(1,length(x)),Cidx',num2cell(krange),'UniformOutput',false));
jvec(1:nz_new) = cell2mat(Cidx');
svec(1:nz_new) = cell2mat(Cdist');
nnz = nz_new;

progress = progress + progress_tick;
waitbar(progress,h_waitbar,sprintf('Computing distance matrix, nnz=%d',nnz));
krange = unique(min(krange + kstep,N));
while progress<1
    [Cidx,Cdist] = rangesearch(x,x(krange,:),1/ep);
    nz_new = sum(cellfun(@length,Cidx));
    ivec(nnz+1:nnz+nz_new) = cell2mat(cellfun(@(x,k)k*ones(1,length(x)),Cidx',num2cell(krange),'UniformOutput',false));
    jvec(nnz+1:nnz+nz_new) = cell2mat(Cidx');
    svec(nnz+1:nnz+nz_new) = cell2mat(Cdist');
    nnz = nnz + nz_new;
    
    krange = unique(min(krange + kstep,N));
    progress = progress + progress_tick;
    waitbar(progress,h_waitbar,sprintf('Computing distance matrix, nnz=%d',nnz));
end
% Add something nonzero to avoid the squeeze Matlab would otherwise apply
svec = svec + eps;
DM = sparse(ivec(1:nnz),jvec(1:nnz),svec(1:nnz),N,N,nnz);
close(h_waitbar);

% Conduct a regression on the data
% This is a vector valued function [r g b], so we need to do each
% dimension separately, although they could be considered together
Xr = zeros(size(Xs));
for k=1:3
    
end