% Density_plot_heatmap

% points=randi([1 10],20,20)
% points=rand(20,20);
clear all;
close all;

path='I:\optogenetics\';
pathin=[path,'Histology\'];
% mousenames=[1 5 6 7 8 9 15 16 17];
mousenames=[1 5 6 7 8 9 12 15 16 17 18 19 20 21];
numanim=length(mousenames);


for nn = 1:numanim
ftips=[-1.15 0.25 -5.2;... %1
       -0.10 0.80 -5.2;... %5
       -0.45 0.70 -5.2;... %6
       -0.90 0.70 -4.9;... %7
       -0.40 0.40 -5.0;... %8
       -0.95 1.05 -5.3;... %9
       -0.45 0.40 -5.6;... %12
        0.10 0.45 -5.7;... %15
        0.00 0.65 -5.4;... %16
        0.10 0.75 -5.6;... %17
        0.00 0.35 -5.6;... %18
        0.10 0.55 -5.5;... %19
        0.25 0.55 -5.5;... %20
        0.10 0.55 -5.4];   %21

%%%% DV - bottom:0
ftips=ftips+[0 0 5.8];

% ftips1=ftips(13,:); %20
ftips2=ftips([2 8:14],:);
% ftips2=ftips([2 8:12 14],:);
ftips3=ftips(5,:); %8
ftips4=ftips(7,:); %12
ftips5=ftips(3,:); %6
ftips6=ftips(6,:); %9
ftips7=ftips(4,:); %7
ftips8=ftips(1,:); %1

%%%%%%%% AP axis reversed

end

% x1=ftips1(:,1); y1=ftips1(:,2); z1=ftips1(:,3);
x2=ftips2(:,1); y2=ftips2(:,2); z2=ftips2(:,3);
x3=ftips3(:,1); y3=ftips3(:,2); z3=ftips3(:,3);
x4=ftips4(:,1); y4=ftips4(:,2); z4=ftips4(:,3);
x5=ftips5(:,1); y5=ftips5(:,2); z5=ftips5(:,3);
x6=ftips6(:,1); y6=ftips6(:,2); z6=ftips6(:,3);
x7=ftips7(:,1); y7=ftips7(:,2); z7=ftips7(:,3);
x8=ftips8(:,1); y8=ftips8(:,2); z8=ftips8(:,3);




% sky blue: 86/255 180/255 233/255 (LPO)
% bluish green: 0 158/255 115/255 (nLPO)

% blue: 0 114/255 178/255 (Wake)
% orange: 230/255 159/255 0 (NREM)
% Reddish purple: 204/255 121/255 167/255 (sedation)
% Vermillion: 213/255 94/255 0 (REM)

% sky blue 25% transparent: 86/255 180/255 233/255 (2-min photostim)


% ds=100; % dot size
cs=2700; % circle size
figure
% stem3(x1, y1, z1,'filled');
stem3(x2, y2, z2, 'filled','Color',[86/255 180/255 233/255]);
hold on
stem3(x3, y3, z3, 'filled','Color',[0 158/255 115/255]);
stem3(x4, y4, z4, 'filled','Color',[0 158/255 115/255]);
stem3(x5, y5, z5, 'filled','Color',[0 158/255 115/255]);
stem3(x6, y6, z6, 'filled','Color',[0 158/255 115/255]);
stem3(x7, y7, z7, 'filled','Color',[0 158/255 115/255]);
stem3(x8, y8, z8, 'filled','Color',[0 158/255 115/255]);
hold on


ax = gca;
ax.ColorOrderIndex = 1;
ax.FontSize = 10;
% h1=scatter3(x1, y1, z1, cs);
% hold on
h2=scatter3(x2, y2, z2, cs,[86/255 180/255 233/255]);
hold on
h3=scatter3(x3, y3, z3, cs,[0 158/255 115/255]);
h4=scatter3(x4, y4, z4, cs,[0 158/255 115/255]);
h5=scatter3(x5, y5, z5, cs,[0 158/255 115/255]);
h6=scatter3(x6, y6, z6, cs,[0 158/255 115/255]);
h7=scatter3(x7, y7, z7, cs,[0 158/255 115/255]);
h8=scatter3(x8, y8, z8, cs,[0 158/255 115/255]);

% hold on
% stem3(0, 0.62, 0, '*r'); % VLPO ML:0.62
% h0=scatter3(0, 0.62, 0, 6400,'r');  % VLPO ML:0.62

view(187,27)
set (gca,'color','none')
% set(gca,'XTick',[],'FontSize',12)%, xticklabels({})
% set(gca,'YTick',[],'FontSize',12)%, yticklabels({})
% set(gca,'ZTick',[],'FontSize',12)%, zticklabels({})
xlim([-1.4 0.601])
ylim([0.0 1.2])
zlim([0 1])
ax.TickLength = [0.02 0.016];
ax.XTick = [-1.4:0.2:0.4];
ax.XTickLabel = {'','-1.2','','-0.8','','-0.4','','0','','0.4',''}
ax.YTick = [0.0:0.2:1.2];
ax.YTickLabel = {'0.0','','0.4','','0.8','','1.2'}
ax.ZTick = [0:0.2:1];
ax.ZTickLabel = {'0','0.2','0.4','0.6','0.8','1.0'}
% ax.XMinorTick = 'on'
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16)
xlabel('AP [mm]','FontSize',18,'Rotation',-2.05)
ylabel('ML [mm]','FontSize',18,'Rotation',70)
zlabel('DV from the bottom [mm]','FontSize',18)
x.FontSmoothing = 'on'
% grid on
% axis off

% 	MarkerFaceAlpha
% MarkerEdgeAlpha
% 
% 
% 
% h1=scatter3(x1, y1, z1, ds, 'filled');
% % h1.MarkerFaceColor = [.7 .1 .7];
% hold on
% h2=scatter3(x2, y2, z2, ds,'filled');
% % h2.MarkerFaceColor = [0 0.2 1];
% h3=scatter3(x3, y3, z3, ds, 'filled');
% % h3.MarkerFaceColor = [0 0.7 0.8];
% h4=scatter3(x4, y4, z4, ds, 'filled');
% % h4.MarkerFaceColor = [0.4 1 0.1];
% h5=scatter3(x5, y5, z5, ds, 'filled');
% % h5.MarkerFaceColor = [1 0.2 0.0];
% h6=scatter3(x6, y6, z6, ds, 'filled');
% h7=scatter3(x7, y7, z7, ds, 'filled');
% h8=scatter3(x8, y8, z8, ds, 'filled');



% 
% xticklabels = -0.4:0.2:1.5;
% xticks = linspace(1, size(Ksum, 2), numel(xticklabels));
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
% 
% yticklabels = [-1.5:0.2:0.4];
% yticks = linspace(1, size(Ksum, 1), numel(yticklabels));
% set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))
% colorbar
% colorbar('Ticks',[0, 1, 2, 3, 4],...
%          'TickLabels',{'0','1','2','3','4'})
%      
%      
% figname=['DensityPlot_GDCh_posi_8anim'];
% saveas(gcf,[pathin,figname],'tiff')
% 
% 
% mousenames2=[3 10 11 14];
% numanim=length(mousenames2);
% b = ones(40, 40);
% KsF=[];

% 
% figure
% K3D=reshape(KsF,[200,200,nn]);
% KsumF=sum(K3D,3);
% eval(['save ',pathin,'GDCh_false.dat KsF KsumF -ascii']);
% imagesc(KsumF);
% 
% xticklabels = -0.4:0.1:1.5;
% xticks = linspace(1, size(KsumF, 2), numel(xticklabels));
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
% 
% yticklabels = [-1.5:0.1:0.4];
% yticks = linspace(1, size(KsumF, 1), numel(yticklabels));
% set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))
% colorbar
