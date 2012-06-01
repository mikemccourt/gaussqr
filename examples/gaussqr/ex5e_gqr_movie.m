% ex5e_gqr_movie
function ex5e_gqr_movie(sol_store,t_store,N,x)
X = reshape(x(:,1),N,N);
Y = reshape(x(:,2),N,N);

window_size = [122 204 971 504];
movie_figure = figure('Position',window_size);

numframes = length(sol_store);
%     movie_matrix = moviein(numframes,movie_figure,window_size);
set(movie_figure,'NextPlot','replacechildren')
for k=1:numframes
    u = sol_store{k};
    A = reshape(u(1:N*N),N,N);
    H = reshape(u(1+N*N:end),N,N);

    subplot(1,2,1)
    surf(X,Y,A,'EdgeColor','none') ,title(sprintf('Activator, t=%d',t_store(k)))
    zlim([0,.1]),shading interp,grid off,xlabel('x'),ylabel('y'),axis square
    subplot(1,2,2)
    surf(X,Y,H,'EdgeColor','none'),title(sprintf('Inhibitor, t=%d',t_store(k)))
    zlim([0,.1]),shading interp,grid off,xlabel('x'),ylabel('y'),axis square
    cmap = colormap('jet');colormap(flipud(cmap))
    pause(.1)
    %         movie_matrix(:,k) = getframe(movie_figure,window_size);
    %         movie_matrix(:,k) = getframe(movie_figure);
end