function plotfcn(output)

figure()
plot3(output.x(1,:),output.x(2,:),output.x(3,:),'LineWidth',1.5);
hold on
grid on
% plot3(output.xref(1,:),output.xref(2,:),output.xref(3,:),'LineWidth',1.5,'LineStyle','--','Color',[1,0.75,0])
quiver3(0,0,0,20,0,0,'LineWidth',1.5,'Color',[1,0,0]);
quiver3(0,0,0,0,20,0,'LineWidth',1.5,'Color',[0,1,0]);
quiver3(0,0,0,0,0,20,'LineWidth',1.5,'Color',[0,0,1]);
phobos = stlread('Phobos_200k.stl');
patch('Vertices',phobos.Vertices,'Faces',phobos.Faces,'FaceColor',[0.455,0.380,0.365],'EdgeAlpha',0);
light
ax = gca();
ax.DataAspectRatio = [1,1,1];
hold off

figure()
plot(output.t,output.rnorm,'LineWidth',1.5)

% figure()
% plot(output.t,output.u,'LineWidth',1.5)
% grid on
% title('Control action')

% figure()
% plot(output.t,output.unorm,'LineWidth',1.5)
% grid on
% title('Control effort')

end
