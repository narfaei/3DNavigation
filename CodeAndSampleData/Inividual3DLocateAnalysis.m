%% Locate 3D Individual streemlined Analysis
clear
clc
close all
%% load Subject data
MoveHead3D= xlsread('Data-1-3DHead.csv');
MoveReach3D= xlsread('Data-2-3DReach.csv');
PinnedReach = xlsread('Data-3-PinnedReach.csv');
Controller3DHead= xlsread('Data-4-Controller3DHead.csv');
set(0, 'defaultFigureRenderer', 'painters');
%% Define Colors

myBlue=[0,0.447000000000000,0.741000000000000];
myCyan=[0.301000000000000,0.745000000000000,0.933000000000000];
myPurple=[0.494000000000000,0.184000000000000,0.556000000000000];
myOrange=[0.850000000000000,0.325000000000000,0.0980000000000000];
myRed=[0.635000000000000,0.0780000000000000,0.184000000000000];
myGreen= [0.466000000000000,0.674000000000000,0.188000000000000];

mycolormap = customcolormap([0 0.7 0.9 1], [myRed; myPurple;myCyan; myBlue; ]);
set(0, 'defaultFigureRenderer', 'painters');
%% cleaning and formating data
PinnedReach(:,101:end)=[];
MoveHead3D(:,101:end)=[];
MoveReach3D(:,101:end)=[];
Controller3DHead(:,101:end)=[];

PinnedReach=PinnedReach';
MoveHead3D=MoveHead3D';
MoveReach3D=MoveReach3D';
Controller3DHead=Controller3DHead';

PinnedReachError=abs(PinnedReach(:,2:4) - PinnedReach(:,5:7));
MoveHead3DError=abs(MoveHead3D(:,2:4) - MoveHead3D(:,5:7));
MoveReach3DError=abs(MoveReach3D(:,2:4) - MoveReach3D(:,5:7));
Controller3DHeadError=abs(Controller3DHead(:,2:4) - Controller3DHead(:,5:7));

CartPinnedReachError=(PinnedReach(:,2:4) - PinnedReach(:,5:7));
CartMoveHead3DError=(MoveHead3D(:,2:4) - MoveHead3D(:,5:7));
CartMoveReach3DError=(MoveReach3D(:,2:4) - MoveReach3D(:,5:7));
CartController3DHeadError=(Controller3DHead(:,2:4) - Controller3DHead(:,5:7));

mean(abs(CartController3DHeadError))
mean(PinnedReachError)
mean(PinnedReachError)
sum(Controller3DHeadError)

%removing outliers based on their error (more than 3std from the mean is removed)
[CleanedPinnedReachError,reachErrorindecies] = rmoutliers(PinnedReachError,'mean');
[CleanedMoveHead3DError,MoveHead3Dindecies] = rmoutliers(MoveHead3DError,'mean');
[CleanedMoveReach3DError,Movereach3Dindecies] = rmoutliers(MoveReach3DError,'mean');
[CleanedController3DHeadError,Controller3Dindecies] = rmoutliers(Controller3DHeadError,'mean');

mean(abs(CleanedController3DHeadError))

PinnedReachError= mean(CleanedPinnedReachError);
MoveHead3DError=mean(CleanedMoveHead3DError);
MoveReach3DError=mean(CleanedMoveReach3DError);
Controller3DHeadError=mean(CleanedController3DHeadError);

% from the main data set take out the outliers
PinnedReach(reachErrorindecies,:)=[];
MoveHead3D(MoveHead3Dindecies,:)=[];
MoveReach3D(Movereach3Dindecies,:)=[];
Controller3DHead(Controller3Dindecies,:)=[];

PinnedReachSTE=PinnedReach(:,2:4) - PinnedReach(:,5:7);
MoveHead3DSTE=MoveHead3D(:,2:4) - MoveHead3D(:,5:7);
MoveReach3DSTE=MoveReach3D(:,2:4) - MoveReach3D(:,5:7);
Controller3DHeadSTE=Controller3DHead(:,2:4) - Controller3DHead(:,5:7);

PinnedReachSTE= std(PinnedReachSTE) / sqrt(length(PinnedReachSTE));
MoveHead3DSTE= std(MoveHead3DSTE) / sqrt(length(MoveHead3DSTE));
MoveReach3DSTE= std(MoveReach3DSTE) / sqrt(length(MoveReach3DSTE));
Controller3DHeadSTE= std(Controller3DHeadSTE) / sqrt(length(Controller3DHeadSTE));

PinnedReachSTE = PinnedReachSTE(:,[1 3 2]);
MoveHead3DSTE = MoveHead3DSTE(:,[1 3 2]);
MoveReach3DSTE = MoveReach3DSTE(:,[1 3 2]);
Controller3DHeadSTE = Controller3DHeadSTE(:,[1 3 2]);

x=[1,2,3];

[h,p,ci,stats] =ttest2(CleanedController3DHeadError(:,1),CleanedController3DHeadError(:,2));
[h,p,ci,stats] =ttest2(CleanedController3DHeadError(:,3),CleanedController3DHeadError(:,2));


%% bar plot error in dimensions with error bars
figure(1)

subplot(1,4,1)
bar(x,[PinnedReachError(1),PinnedReachError(3),PinnedReachError(2)]*100)     
xticklabels({'X','Y','Z'})
hold on
errorbar(x,[PinnedReachError(1),PinnedReachError(3),PinnedReachError(2)]*100,PinnedReachSTE*100,'.','Color','k'); 
ylim([0,170])
grid on
title('Reach')
xlabel('Dimension')
ylabel('Error in cm')

subplot(1,4,2)
bar(x,[MoveReach3DError(1),MoveReach3DError(3),MoveReach3DError(2)]*100)     
xticklabels({'X','Y','Z'})
hold on
errorbar(x,[MoveReach3DError(1),MoveReach3DError(3),MoveReach3DError(2)]*100, MoveReach3DSTE*100,'.','Color','k'); 
ylim([0,170])
grid on
title('Hybrid Reach')
xlabel('Dimension')
ylabel('Error in cm')

subplot(1,4,3)
bar(x,[MoveHead3DError(1),MoveHead3DError(3),MoveHead3DError(2)]*100)     
xticklabels({'X','Y','Z'})
hold on
errorbar(x,[MoveHead3DError(1),MoveHead3DError(3),MoveHead3DError(2)]*100,MoveHead3DSTE*100,'.','Color','k'); 
ylim([0,170])
grid on
title('Hybrid Navigation')
xlabel('Dimension')
ylabel('Error in cm')


subplot(1,4,4)
bar(x,[Controller3DHeadError(1),Controller3DHeadError(3),Controller3DHeadError(2)]*100)     
xticklabels({'X','Y','Z'})
hold on
errorbar(x,[Controller3DHeadError(1),Controller3DHeadError(3),Controller3DHeadError(2)]*100,Controller3DHeadSTE*100,'.','Color','k'); 
ylim([0,170])
grid on
title('Virtual Navigation')
xlabel('Dimension')
ylabel('Error in cm')

s = pwd;
name=string(s)+'\SummaryStats\MainErrorBarPlot.jpg';
saveas(gcf,name)

% run pairwise ttest to see if difference is significant
[h,p,ci,stats] =ttest2(CleanedController3DHeadError(:,1),CleanedController3DHeadError(:,2));
[h,p,ci,stats] =ttest2(CleanedController3DHeadError(:,3),CleanedController3DHeadError(:,2));
% run pairwise ttest to see if x and y are the same
[h,p,ci,stats] =ttest2(CleanedController3DHeadError(:,1),CleanedController3DHeadError(:,3));

[h,p,ci,stats] =ttest2(CleanedMoveHead3DError(:,1),CleanedMoveHead3DError(:,2));
[h,p,ci,stats] =ttest2(CleanedMoveHead3DError(:,3),CleanedMoveHead3DError(:,2));
% run pairwise ttest to see if x and y are the same
[h,p,ci,stats] =ttest2(CleanedMoveHead3DError(:,1),CleanedMoveHead3DError(:,3));

[h,p,ci,stats] =ttest2(CleanedMoveReach3DError(:,1),CleanedMoveReach3DError(:,2));
[h,p,ci,stats] =ttest2(CleanedMoveReach3DError(:,3),CleanedMoveReach3DError(:,2));
% run pairwise ttest to see if x and y are the same
[h,p,ci,stats] =ttest2(CleanedMoveReach3DError(:,1),CleanedMoveReach3DError(:,3));

[h,p,ci,stats] =ttest2(CleanedPinnedReachError(:,1),CleanedPinnedReachError(:,2));
[h,p,ci,stats] =ttest2(CleanedPinnedReachError(:,3),CleanedPinnedReachError(:,2));
% run pairwise ttest to see if x and y are the same
[h,p,ci,stats] =ttest2(CleanedPinnedReachError(:,1),CleanedPinnedReachError(:,3));

% Save summary and plots
ErrorSummaryData(1,:)=PinnedReachError;
ErrorSummaryData( 1, end+1:end+3 )=MoveReach3DError;
ErrorSummaryData( 1, end+1:end+3 )=MoveHead3DError;
ErrorSummaryData( 1, end+1:end+3 )=Controller3DHeadError;
ErrorSummaryData( 1, end+1:end+3 )=PinnedReachSTE;
ErrorSummaryData( 1, end+1:end+3 )=MoveReach3DSTE;
ErrorSummaryData( 1, end+1:end+3 )=MoveHead3DSTE;
ErrorSummaryData( 1, end+1:end+3 )=Controller3DHeadSTE;

colNames = {'PinnedMeanErrorX','PinnedMeanErrorZ','PinnedMeanErrorY','MoveReachErrorX','MoveReachErrorZ','MoveReachErrorY','MoveHeadErrorX','MoveHeadErrorY','MoveHeadErrorZ','ControllerErrorX','ControllerErrorZ','ControllerErrorY','PinnedSTEX','PinnedSTEZ','PinnedSTEY','MoveReachSTEX','MoveReachSTEZ','MoveReachSTEY','MoveHeadSTEX','MoveHeadSTEY','MoveHeadSTEZ','ControllerSTEX','ControllerSTEZ','ControllerSTEY'};
rowNames = {'1'};
sTable = array2table(ErrorSummaryData,'RowNames',rowNames,'VariableNames',colNames);

name=string(s)+'\SummaryStats\ErrorSummaryData.csv';
% Save it
writetable(sTable,name);

% save the errors in seperate excel
CartError(:,1)=CartController3DHeadError(:,1);
CartError(:,2)=CartController3DHeadError(:,3);
CartError(:,3)=CartController3DHeadError(:,2);

sTable = array2table(CartError');
name=string(s)+'\SummaryStats\CartesianError.csv';
writetable(sTable,name);
name=string(s)+'\SummaryStats\CartesianError.svg';
writetable(sTable,name);


%% Violin Plot

x0=400;
y0=200;
width=200;
height=600;
wWidth=1.8;
LfontSize= 12;
TfontSize= 16;
ticksizeX= 12;
ticksizeY= 12;
DotSize=150;
MeanSize= 25;
SideLineW=1;
colorBlue=[.4,.4,.8];
x=[.5,1.5,2.5];
Csize= 11;
 
figure(2)
subplot(1,4,1)
location=ones(3,length(CleanedPinnedReachError));
location=x.*location';
violins= violinplot([CleanedPinnedReachError(:,1),CleanedPinnedReachError(:,3),CleanedPinnedReachError(:,2)]*100,'ShowMean' ,'MedianColor',[0,0,0],'BoxWidth' ,0.0,'Width',.2);
xticklabels({'X','Y','Z'});
errorbar([PinnedReachError(1),PinnedReachError(3),PinnedReachError(2)]*100, PinnedReachSTE*100,'.','Color','k','MarkerSize',MeanSize,'LineWidth',wWidth,'CapSize',Csize); 
ylim([0,150]);
xlim([.5,3.5]);
title('Reach','fontweight','bold','FontSize',TfontSize);
xlabel('Dimension','fontweight','bold','FontSize',LfontSize);
ylabel('Error in cm','fontweight','bold','FontSize',LfontSize);
set(gca, 'YGrid', 'off', 'XGrid', 'on');
ax = gca;
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
ax.LineWidth = SideLineW;

subplot(1,4,2)
location=ones(3,length(CleanedMoveReach3DError));
location=x.*location';

violins= violinplot([CleanedMoveReach3DError(:,1),CleanedMoveReach3DError(:,3),CleanedMoveReach3DError(:,2)]*100,'ShowMean' ,'MedianColor',[0,0,0],'BoxWidth' ,0.0,'Width',.2);
xticklabels({'X','Y','Z'});
errorbar([MoveReach3DError(1),MoveReach3DError(3),MoveReach3DError(2)]*100, MoveReach3DSTE*100,'.','Color','k','MarkerSize',MeanSize,'LineWidth',wWidth,'CapSize',Csize); 
ylim([0,150]);
xlim([.5,3.5]);
title('Reach Hybrid','fontweight','bold','FontSize',TfontSize);
xlabel('Dimension','fontweight','bold','FontSize',LfontSize);
ylabel('Error in cm','fontweight','bold','FontSize',LfontSize);
set(gca, 'YGrid', 'off', 'XGrid', 'on')
ax = gca;
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
ax.LineWidth = SideLineW;

subplot(1,4,3)
location=ones(3,length(CleanedMoveHead3DError));
location=x.*location';
violins= violinplot([CleanedMoveHead3DError(:,1),CleanedMoveHead3DError(:,3),CleanedMoveHead3DError(:,2)]*100,'ShowMean' ,'MedianColor',[0,0,0],'BoxWidth' ,0.0,'Width',.2);
xticklabels({'X','Y','Z'});
errorbar([MoveHead3DError(1),MoveHead3DError(3),MoveHead3DError(2)]*100, MoveHead3DSTE*100,'.','Color','k','MarkerSize',MeanSize,'LineWidth',wWidth,'CapSize',Csize); 
ylim([0,150]);
xlim([.5,3.5]);
title('Hybrid','fontweight','bold','FontSize',TfontSize);
xlabel('Dimension','fontweight','bold','FontSize',LfontSize);
ylabel('Error in cm','fontweight','bold','FontSize',LfontSize);
set(gca, 'YGrid', 'off', 'XGrid', 'on');
ax = gca;
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
ax.LineWidth = SideLineW;

subplot(1,4,4)
location=ones(3,length(CleanedController3DHeadError));
location=x.*location';
violins= violinplot([CleanedController3DHeadError(:,1),CleanedController3DHeadError(:,3),CleanedController3DHeadError(:,2)]*100,'ShowMean' ,'MedianColor',[0,0,0],'BoxWidth' ,0.0,'Width',.2);

xticklabels({'X','Y','Z'});
errorbar([Controller3DHeadError(1),Controller3DHeadError(3),Controller3DHeadError(2)]*100, Controller3DHeadSTE*100,'.','Color','k','MarkerSize',MeanSize,'LineWidth',wWidth,'CapSize',Csize); 
ylim([0,150]);
xlim([.5,3.5]);
title('Virtual','fontweight','bold','FontSize',TfontSize);
xlabel('Dimension','fontweight','bold','FontSize',LfontSize);
ylabel('Error in cm','fontweight','bold','FontSize',LfontSize);
set(gca, 'YGrid', 'off', 'XGrid', 'on');
ax = gca;
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
ax.LineWidth = SideLineW;

x0=400;
y0=200;
width=1000;
height=400;
set(gcf,'position',[x0,y0,width,height]);


s = pwd;
name=string(s)+'\SummaryStats\CartesianErrorViolinPlot.jpg';
saveas(gcf,name)
name=string(s)+'\SummaryStats\CartesianErrorViolinPlot.svg';
saveas(gcf,name)


%% Error Fields: Plot errors compared to the actual location of stimuli (stim location as [0,0,0])

wWidth=2.5;
LfontSize= 22;
TfontSize= 24;
ticksizeX= 14;
ticksizeY= 14;
ticksizeZ= 14;
DotSize=40;

figure(3)
scatter3(Controller3DHead(:,5)-Controller3DHead(:,2),Controller3DHead(:,7)-Controller3DHead(:,4),Controller3DHead(:,6)-Controller3DHead(:,3),'filled','SizeData',DotSize)
hold on

xdiffCentered=zeros(length(Controller3DHead),1);
xdiffCentered(:,2)=Controller3DHead(:,5)-Controller3DHead(:,2);
ydiffCentered=zeros(length(Controller3DHead),1);
ydiffCentered(:,2)=Controller3DHead(:,7)-Controller3DHead(:,4);
zdiffCentered=zeros(length(Controller3DHead),1);
zdiffCentered(:,2)=Controller3DHead(:,6)-Controller3DHead(:,3);
xdiffCentered=xdiffCentered';
ydiffCentered=ydiffCentered';
zdiffCentered=zdiffCentered';

for ii=1:length(zdiffCentered)
plot3(xdiffCentered(:,ii),ydiffCentered(:,ii),zdiffCentered(:,ii), 'linewidth',.2,'color','k')
quiver3(xdiffCentered(1,ii),ydiffCentered(1,ii),zdiffCentered(1,ii),xdiffCentered(2,ii),ydiffCentered(2,ii),zdiffCentered(2,ii), 'linewidth',.2,'color','k')
end
axis equal
title('3D Error from target')
xlabel('x')
ylabel('y')
zlabel('z')

axis equal
title('Error from target for one participant','FontSize',TfontSize)
xlabel('x (m)','FontSize',LfontSize)
ylabel('y (m)','FontSize',LfontSize)
zlabel('z (m)','FontSize',LfontSize)
xticks([-4 , -3, -2 , -1, 0 , 1, 2, 3, 4])
yticks([-4 , -3, -2 , -1, 0 , 1, 2, 3, 4])
zticks([-4 , -3, -2 , -1, 0 , 1, 2, 3, 4])
xlim([-4 4])
ylim([-4 4])
zlim([-4 4])
ax = gca;
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
ax.ZAxis.FontSize =ticksizeZ;
axis equal
x0=10;
y0=10;
width=500;
height=500;
set(gcf,'position',[x0,y0,width,height])
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
set(get(gca, 'ZAxis'), 'FontWeight', 'bold');
ax.LineWidth = SideLineW;

xlim([-1 1])
ylim([-1 1])
zlim([-2 2])
xticks([-4 , -3, -2 , -1, 0 , 1, 2, 3, 4])
yticks([-4 , -3, -2 , -1, 0 , 1, 2, 3, 4])
zticks([-4 , -3, -2 , -1, 0 , 1, 2, 3, 4])

x0=10;
y0=10;
width=500;
height=500;
set(gcf,'position',[x0,y0,width,height])
name=string(s)+'\SummaryStats\ErrorField.jpg';
saveas(gcf,name)
name=string(s)+'\SummaryStats\ErrorField.svg';
saveas(gcf,name)

%% Calculate central tendency and dispersion

figure(4)
diffCentered=cat(1, (xdiffCentered(2,:)) , (ydiffCentered(2,:)) ,(zdiffCentered(2,:))) ;
CT=mean(diffCentered,2);
CovM=cov(diffCentered');
imagesc(CovM);

set(gca, 'XTick', 1:3); % center x-axis ticks on bins
set(gca, 'YTick', 1:3); % center y-axis ticks on bins
set(gca, 'XTickLabel', {'X' ,'Y' ,'Z'}); % set x-axis labels
set(gca, 'YTickLabel',  {'X' ,'Y' ,'Z'}); % set y-axis labels
title('Covariance of Error from veridical', 'FontSize', 10); % set title
colormap(jet(256));% Choose jet or any other color scheme
colorbar;

s = pwd;
name=string(s)+'\SummaryStats\CovariancePlot.jpg';
saveas(gcf,name)

% Save Data
CovData(1,1)=CovM(1,1);
CovData(1,2)=CovM(2,2);
CovData(1,3)=CovM(3,3);
CovData(1,4)=CovM(1,2);
CovData(1,5)=CovM(1,3);
CovData(1,6)=CovM(2,3);
CovData(1,7:9)=CT';

colNames = {'VarX','VarY','VarZ','CovXY','CovXZ','CovYZ','meanX','meanY','meanZ'};
rowNames = {'1'};
sTable = array2table(CovData,'RowNames',rowNames,'VariableNames',colNames);
name=string(s)+'\SummaryStats\CovarianceData.csv';
writetable(sTable,name);
name=string(s)+'\SummaryStats\CovarianceData.svg';
writetable(sTable,name);

%% Plot error vectors in 3D controller
wWidth=1.8;
LfontSize= 20;
TfontSize= 20;
ticksizeX= 14;
ticksizeY= 14;
DotSize=35;
MeanSize= 25;
SideLineW=1;
myAlpha=0.9;

figure(5)
subplot(1,3,1)
s=scatter(Controller3DHead(:,5), Controller3DHead(:,2),'filled','SizeData',DotSize,'MarkerFaceColor',myBlue,'MarkerFaceAlpha',myAlpha,'MarkerEdgeAlpha',myAlpha);
ax = gca;
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
ax.LineWidth = SideLineW;
title('X','fontweight','bold','FontSize',TfontSize)
xlabel('Navigated(m)','fontweight','bold','FontSize',LfontSize)
ylabel('Actual(m)','fontweight','bold','FontSize',LfontSize)
hold on
plot([-4 4],[-4 4],'--','LineWidth',2,'color','k')
axis equal
xlim([-2 2])
ylim([-2 2])
p = polyfit(Controller3DHead(:,5), Controller3DHead(:,2),1);
x1 = linspace(-4,4,2);
y1 = polyval(p,x1);
plot(x1,y1,'LineWidth',2,'color','r')
ax.LineWidth = SideLineW;
xticks([-4 , -3, -2 , -1, 0 , 1, 2, 3, 4])
yticks([-4 , -3, -2 , -1, 0 , 1, 2, 3, 4])
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
slopex=abs(y1(1)/y1(2));
interceptx=[ -p(2) / p(1),p(2)];

subplot(1,3,2)
scatter(Controller3DHead(:,7), Controller3DHead(:,4),'filled','SizeData',DotSize,'MarkerFaceColor',myCyan,'MarkerFaceAlpha',myAlpha,'MarkerEdgeAlpha',myAlpha)
axis equal
ax = gca;
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
ax.LineWidth = SideLineW;
hold on
plot([-4 4],[-4 4],'--','LineWidth',2,'color','k')
axis equal
xlim([-2 2])
ylim([-2 2])
p = polyfit(Controller3DHead(:,7), Controller3DHead(:,4),1);
x1 = linspace(-4,4,2);
y1 = polyval(p,x1);
plot(x1,y1,'LineWidth',2,'color','r')
ax.LineWidth = SideLineW;
title('Y','fontweight','bold','FontSize',TfontSize)
xlabel('Navigated(m)','fontweight','bold','FontSize',LfontSize)
ylabel('Actual(m)','fontweight','bold','FontSize',LfontSize)
xticks([-4 , -3, -2 , -1, 0 , 1, 2, 3, 4])
yticks([-4 , -3, -2 , -1, 0 , 1, 2, 3, 4])
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
width=1000;
height=300;
x0=100;
y0=100;
set(gcf,'position',[x0,y0,width,height])
slopey=abs(y1(1)/y1(2));
intercepty=[ -p(2) / p(1),p(2)];

subplot(1,3,3)
scatter(Controller3DHead(:,6), Controller3DHead(:,3),'filled','SizeData',DotSize,'MarkerFaceColor',myPurple,'MarkerFaceAlpha',myAlpha,'MarkerEdgeAlpha',myAlpha)
ax = gca;
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
ax.LineWidth = SideLineW;
hold on
plot([-4 7],[-4 7],'--','LineWidth',2,'color','k')
axis equal
xlim([2 6])
ylim([2 6])
p = polyfit(Controller3DHead(:,6), Controller3DHead(:,3),1);
x1 = linspace(-4,8,2);
y1 = polyval(p,x1);
plot(x1,y1,'LineWidth',2,'color','r')
ax.LineWidth = SideLineW;
title('Z','fontweight','bold','FontSize',TfontSize)
xlabel('Navigated (m)','fontweight','bold','FontSize',LfontSize)
ylabel('Actual(m)','fontweight','bold','FontSize',LfontSize)
xticks([-4 , -3, -2 , -1, 0 , 1, 2, 3, 4 ,5,6])
yticks([-4 , -3, -2 , -1, 0 , 1, 2, 3, 4,5,6])
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
slopez=abs(y1(1)/y1(2));
interceptz=[ -p(2) / p(1),p(2)];


s = pwd;
name=string(s)+'\SummaryStats\CorrellationPlot.jpg';
saveas(gcf,name)
name=string(s)+'\SummaryStats\CorrellationPlot.svg';
saveas(gcf,name)

%%
figure(6)

scatter3(Controller3DHead(:,2),Controller3DHead(:,4),Controller3DHead(:,3), 'filled')
hold on
scatter3(Controller3DHead(:,5),Controller3DHead(:,7),Controller3DHead(:,6),'filled')

xdiff=Controller3DHead(:,2);
xdiff(:,2)=Controller3DHead(:,5);
ydiff=Controller3DHead(:,4);
ydiff(:,2)=Controller3DHead(:,7);
zdiff=Controller3DHead(:,3);
zdiff(:,2)=Controller3DHead(:,6);
xdiff=xdiff';
ydiff=ydiff';
zdiff=zdiff';
title("Trajectory of error compared to location of stimuli")
for ii=1:length(xdiff)
plot3(xdiff(:,ii),ydiff(:,ii),zdiff(:,ii), 'linewidth',1.5,'color','k')
quiver3(xdiff(1,ii),ydiff(1,ii),zdiff(1,ii),xdiff(2,ii)-xdiff(1,ii),ydiff(2,ii)-ydiff(1,ii),zdiff(2,ii)-zdiff(1,ii), 'linewidth',1.5,'color','k','MaxHeadSize',5)
end

xlabel("X")
ylabel("Y")
zlabel("Z")
width=800;
height=1000;
x0=100;
y0=100;
set(gcf,'position',[x0,y0,width,height])

name=string(s)+'\SummaryStats\ErrorInSpacePlot.jpg';
saveas(gcf,name)
name=string(s)+'\SummaryStats\ErrorInSpacePlot.svg';
saveas(gcf,name)

%% Timecourse of the errors

xerror= Controller3DHead(:,5)-Controller3DHead(:,2);
yerror=Controller3DHead(:,7)-Controller3DHead(:,4);
zerror=Controller3DHead(:,6)-Controller3DHead(:,3);

trialNum=Controller3DHead(:,1);

figure(7)
subplot(3,1,1)
plot(trialNum,xerror)
title("Error as a function of trial number")
ylim([-3,3])
yline(0)
grid on
ylabel('xError')
subplot(3,1,2)
plot(trialNum,yerror)
ylim([-3,3])
grid on
ylabel('yError')
yline(0)
subplot(3,1,3)
plot(trialNum,zerror)
ylim([-3,3])
grid on
ylabel('zError')
yline(0)

name=string(s)+'\SummaryStats\ErrorInTimePlot.jpg';
saveas(gcf,name)
name=string(s)+'\SummaryStats\ErrorInTimePlot.svg';
saveas(gcf,name)
%% Skewness of data
CartErrorXall=CartController3DHeadError(:,1);
CartErrorYall=CartController3DHeadError(:,3);
CartErrorZall=CartController3DHeadError(:,2);
Step=.05;

figure(8)
subplot(3,1,1)
NumBins=floor((max(CartErrorXall)-min(CartErrorXall))/Step);
histfit(CartErrorXall,NumBins,'kernel')
title('X dimension')
xlabel('error from veridical(m)')
ylabel('Frequency of error')
grid on
xlim([-4,4])
ylim([0,15])
subplot(3,1,2)
NumBins=floor((max(CartErrorYall)-min(CartErrorYall))/Step);
histfit(CartErrorYall,NumBins,'kernel');
xlim([-4,4])
ylim([0,15])
title('Y dimension')
xlabel('error from veridical(m)')
ylabel('Frequency of error')
grid on
subplot(3,1,3)
NumBins=floor((max(CartErrorZall)-min(CartErrorZall))/Step);
histfit(CartErrorZall,NumBins,'kernel');
xlim([-4,4])
ylim([0,15])
title('Z dimension')
xlabel('error from veridical(m)')
ylabel('Frequency of error')
grid on
s = pwd;
name=string(s)+'\SummaryStats\SkewnessPlot.jpg';
saveas(gcf,name)
%One measure of skewness, called Pearson's first coefficient of skewness, is to subtract the mean from the mode, and then divide this difference by the standard deviation of the data. The reason for dividing the difference is so that we have a dimensionless quantity.
% subtract the mean from the mode, and then divide this difference by the standard deviation of the data.
%When you set flag to 0, skewness corrects for the systematic bias, and the following equation applies:
PCSX=skewness(CartErrorXall,0);
PCSY=skewness(CartErrorYall,0);
PCSZ=skewness(CartErrorZall,0);

PCS(:,1)=PCSX;
PCS(:,2)=PCSY;
PCS(:,3)=PCSZ;

colNames = {'PCSX','PCSY','PCSZ'};
rowNames = {'1'};

sTable = array2table(PCS,'RowNames',rowNames,'VariableNames',colNames);
name=string(s)+'\SummaryStats\Skewness.csv';
writetable(sTable,name);
name=string(s)+'\SummaryStats\Skewness.svg';
writetable(sTable,name);


