function SurfacePlot( evalpnts, faces, fun, titl)
% figure('Color',[1 1 1]);
patch('Vertices',evalpnts,'Faces',faces,'FaceVertexCData',fun,...
      'FaceColor','interp','EdgeColor','none');
colorbar('FontSize',12);
axis('equal'); axis('square'); view(-30,30);
xlabel('x','FontWeight','bold','FontSize',14);
ylabel('y','FontWeight','bold','FontSize',14);
zlabel('z','FontWeight','bold','FontSize',14);
title(titl,'FontWeight','bold','FontSize',14);
grid on;
end

