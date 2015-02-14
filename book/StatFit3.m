% StatFit3.m
% This uses the carsmall dataset in Matlab to conduct some statistical data
% fitting on the impact of:
%     [Acceleration Displacement Horsepower Weight]
% on a car's MPG.

% Load, clean and scale the data
load carsmall
xdirty = [Acceleration Displacement Horsepower Weight];
xstr = {'Acceleration','Displacement','Horsepower','Weight'};
ydirty = MPG;
[x,y,shift,scale] = rescale_data(xdirty,ydirty);
x_mean = mean(x);

% Create the surrogate model using a pre-chosen shape parameter
ep = 4;
rbf = @(e,r) exp(-(e*r).^2);
K = rbf(ep,DistanceMatrix(x,x));
coef = K\y;
SM_eval = @(xx) rbf(ep,DistanceMatrix(xx,x))*coef;

% Evaluate cross-sections of the surrogate model
% Two of the inputs will be fixed to their means
% The other two will be allowed to vary over their input range
NN = 50;
x2d = pick2Dpoints([-1 -1],[1 1],NN*[1;1]);
X1 = reshape(x2d(:,1),NN,NN);
X2 = reshape(x2d(:,2),NN,NN);
mean_vecs = ones(NN^2,1)*x_mean;
subplot_cs = {[1 2],[1 3],[1 4],[2 3],[2 4],[3 4]};
h_cs = figure;
for k=1:length(subplot_cs)
    subplot(3,2,k)
    this_subplot = subplot_cs{k};
    xeval = mean_vecs;
    xeval(:,this_subplot) = x2d;
    yeval = SM_eval(xeval);
    surf(X1,X2,reshape(yeval,NN,NN),'edgecolor','none');
    xlabel(xstr{this_subplot(1)});
    ylabel(xstr{this_subplot(2)});
    zlabel('MPG');
    zlim([0 30])
end

xstr = {'Acc','Disp','HP','Weight'};
% Create some parallel coordinate plots
N = 15;
x1 = pickpoints(-1,1,N);
[X1,X2,X3,X4] = ndgrid(x1);
x4d = [X1(:),X2(:),X3(:),X4(:)];
yfull = SM_eval(x4d);
subplot_pc = {[9 14],[14 21],[21 28],[28 55]};
h_pc = figure;
for k=1:length(subplot_pc)
    subplot(2,2,k);
    this_subplot = subplot_pc{k};
    this_ind = yfull>this_subplot(1) & yfull<this_subplot(2);
    xgood = x4d(this_ind,:);
    groups = ceil(3*(max(xgood(:,1)+1,eps))./2);
    parallelcoords(xgood,'Group',groups,'Label',xstr);
    title(sprintf('MPG between %d and %d',this_subplot))
end



% Consider integrating to a constant 1 if needed (need to work on this
% dblquad(@(x,y)pdf2d_eval([repmat(x',length(y'),1),repmat(y',length(x'),1)]),-1,1,-1,1,1e-8)

pause

% Consider a range of epsilon values over which we study the maximum and
% minimum values that the surrogate model takes in the region
% epvec = logspace(0,1,20);
% maxvec = zeros(size(epvec));
% minvec = zeros(size(epvec));
% opt_opts.Display = 'off';
% k = 1;
% for ep=epvec
%     K = rbf(ep,DistanceMatrix(x,x));
%     coef = K\y;
%     SM_eval = @(xx) rbf(ep,DistanceMatrix(xx,x))*coef;
%     yfull = SM_eval(x4d);
%     [c,minind] = min(yfull)
%     [xminmpg,minvec(k)] = fmincon(@(x)SM_eval(x'),x4d(minind,:)',zeros(4),zeros(4,1),[],[],-ones(4,1),ones(4,1),[],opt_opts);
%     [c,maxind] = max(yfull)
%     [xmaxmpg,maxvec(k)] = fmincon(@(x)-SM_eval(x'),x4d(maxind,:)',zeros(4),zeros(4,1),[],[],-ones(4,1),ones(4,1),[],opt_opts);
%     k = k + 1;
% end
% maxvec = -maxvec;
% h_mm = figure;
% semilogx(epvec,[minvec;maxvec])

% Create a sampling strategy for plotting results in 4D
% Just a tensor product uniform for now
% N = 15;
% x1 = pickpoints(-1,1,N);
% x4v = zeros(prod(size(x4)),4);
% for i=1:N
%     for j=1:N
%         for k=1:N
%             for l=1:N
%                 x4v(i + N*(j-1) + N^2*(k-1) + N^3*(l-1),:) = [x1(i),x1(j),x1(k),x1(l)];
%             end
%         end
%     end
% end
% 
% % Consider the results for fixed Acceleration/Displacement indices
% xx_AD = unique(x4v(:,1:2),'rows');
% DDD = DistanceMatrix(x4v(:,1:2),xx_AD)==0;
% xx_AD_ind = cell(1,N^2);
% for k=1:N^2
%     xx_AD_ind{k} = find(DDD(:,k));
% end
% 
% % Evaluate the full 4D interpolant and fix the AD dimensions
% ep = 4;
% K = rbf(ep,DistanceMatrix(x,x));
% yp = rbf(ep,DistanceMatrix(x4v,x))*(K\y);
% yy_AD = zeros(N^2,N^2);
% H_mat = zeros(N^2,N^2);
% W_mat = zeros(N^2,N^2);
% for k=1:N^2
%     yy_AD(k,:) = yp(xx_AD_ind{k})';
%     H_mat(k,:) = x4v(xx_AD_ind{k},3)';
%     W_mat(k,:) = x4v(xx_AD_ind{k},4)';
% end
% h_AD_GQR = figure;
% surf(reshape(H_mat(1,:),N,N),reshape(W_mat(1,:),N,N),reshape(mean(yy_AD,2),N,N))
% xlabel('Horsepower')
% ylabel('Weight')
% zlabel('MPG')
% title('Averaged over Acceleration/Displacement')
% view([1 1 1])
% pause
% for k=1:N^2
%     surf(reshape(H_mat(k,:),N,N),reshape(W_mat(k,:),N,N),reshape(yy_AD(:,k),N,N))
%     view([1 1 1])
%     zlim([0,40])
%     caxis([0,40])
%     title(sprintf('Acceleration = %g, Displacement = %g',xx_AD(k,1),xx_AD(k,2)))
%     xlabel('Horsepower')
%     ylabel('Weight')
%     zlabel('MPG')
%     pause
% end


% Also, can create a tensor product-ish grid using the AD locations that we
% were given in the data.
% This is for computing the averaging over the input data points.
N = 30;
x1 = pickpoints(-1,1,N);
xAD = x(:,1:2);
NAD = size(xAD,1);
xADv = zeros(NAD*N^2,4);
for i=1:NAD
    for j=1:N
        for k=1:N
            xADv(i + NAD*(j-1) + NAD*N*(k-1),:) = [xAD(i,:),x1(j),x1(k)];
        end
    end
end
ep = 4;
K = rbf(ep,DistanceMatrix(x,x));
yp = rbf(ep,DistanceMatrix(xADv,x))*(K\y);
yy_AD = zeros(N^2,1);
H_vec = zeros(N^2,1);
W_vec = zeros(N^2,1);
for k=1:N^2
    yy_AD(k) = mean(yp(1+(k-1)*NAD:k*NAD));
    H_vec(k) = mean(xADv(1+(k-1)*NAD:k*NAD,3)); % These values are all equal
    W_vec(k) = mean(xADv(1+(k-1)*NAD:k*NAD,4));
end
h_HW_avg = figure;
surf(reshape(H_vec,N,N),reshape(W_vec,N,N),reshape(yy_AD,N,N))
xlabel('Horsepower')
ylabel('Weight')
zlabel('MPG')