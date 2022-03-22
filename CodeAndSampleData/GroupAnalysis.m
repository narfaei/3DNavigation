%% Locate 3D Group streemlined Analysis
clear
clc
close all

%%

myBlue=[0,0.447000000000000,0.741000000000000]
myCyan=[0.301000000000000,0.745000000000000,0.933000000000000]
myPurple=[0.494000000000000,0.184000000000000,0.556000000000000]
myOrange=[0.850000000000000,0.325000000000000,0.0980000000000000]
myRed=[0.635000000000000,0.0780000000000000,0.184000000000000]
myGreen= [0.466000000000000,0.674000000000000,0.188000000000000]

mycolormap = customcolormap([0 0.7 0.9 1], [myRed; myPurple;myCyan; myBlue; ]);

%% load Group data
% Pick the Directory of Combined Data
ErrorsData= xlsread('SummaryErrors.csv');

ErrorsData(:,1)=[];

PinnedReachError=ErrorsData(:,1:3) ;
MoveReach3DError=ErrorsData(:,4:6) ;
MoveHead3DError=ErrorsData(:,7:9) ;
Controller3DHeadError=ErrorsData(:,10:12) ;

PinnedReachSTE= std(PinnedReachError) / sqrt(length(PinnedReachError));
MoveHead3DSTE= std(MoveHead3DError) / sqrt(length(MoveHead3DError));
MoveReach3DSTE= std(MoveReach3DError) / sqrt(length(MoveReach3DError));
Controller3DHeadSTE= std(Controller3DHeadError) / sqrt(length(Controller3DHeadError));

PinnedReachMean=mean(PinnedReachError) ;
MoveReach3DMean=mean(MoveReach3DError);
MoveHead3DMean=mean(MoveHead3DError) ;
Controller3DHeadMean=mean(Controller3DHeadError) ;

PinnedReachSTE = PinnedReachSTE(:,[1 3 2]);
MoveHead3DSTE = MoveHead3DSTE(:,[1 3 2]);
MoveReach3DSTE = MoveReach3DSTE(:,[1 3 2]);
Controller3DHeadSTE = Controller3DHeadSTE(:,[1 3 2]);

x=[1,2,3];


%% bar plot error in dimensions with error bars

figure(1)

subplot(1,4,2)
bar(x,[MoveReach3DMean(1),MoveReach3DMean(3),MoveReach3DMean(2)]*100)     
xticklabels({'X','Y','Z'})
hold on
errorbar(x,[MoveReach3DMean(1),MoveReach3DMean(3),MoveReach3DMean(2)]*100, MoveReach3DSTE*100,'.','Color','k'); 
ylim([0,170])
grid on
title('Reach3D')
xlabel('Dimension')
ylabel('Error in cm')

subplot(1,4,3)
bar(x,[MoveHead3DMean(1),MoveHead3DMean(3),MoveHead3DMean(2)]*100)     
xticklabels({'X','Y','Z'})
hold on
errorbar(x,[MoveHead3DMean(1),MoveHead3DMean(3),MoveHead3DMean(2)]*100,MoveHead3DSTE*100,'.','Color','k'); 
ylim([0,170])
grid on
title('HeadPlacement3D')
xlabel('Dimension')
ylabel('Error in cm')

subplot(1,4,1)
bar(x,[PinnedReachMean(1),PinnedReachMean(3),PinnedReachMean(2)]*100)     
xticklabels({'X','Y','Z'})
hold on
errorbar(x,[PinnedReachMean(1),PinnedReachMean(3),PinnedReachMean(2)]*100,PinnedReachSTE*100,'.','Color','k'); 
ylim([0,170])
grid on
title('ReachPinned')
xlabel('Dimension')
ylabel('Error in cm')

subplot(1,4,4)
bar(x,[Controller3DHeadMean(1),Controller3DHeadMean(3),Controller3DHeadMean(2)]*100)     
xticklabels({'X','Y','Z'})
hold on
errorbar(x,[Controller3DHeadMean(1),Controller3DHeadMean(3),Controller3DHeadMean(2)]*100,Controller3DHeadSTE*100,'.','Color','k'); 
ylim([0,170])
grid on
title('Controller3DHead')
xlabel('Dimension')
ylabel('Error in cm')

s = pwd;
name=string(s)+'\GroupErrorBarPlot.jpg';
saveas(gcf,name)

s = pwd;
name=string(s)+'\GroupErrorBarPlot.svg';
saveas(gcf,name)

%% Stat Analysis
% run pairwise ttests
[h,p,ci,stats] =ttest2(Controller3DHeadError(:,1),Controller3DHeadError(:,2))
[h,p,ci,stats] =ttest2(Controller3DHeadError(:,1),Controller3DHeadError(:,3))

[p,t,stats] = anova1(Controller3DHeadError)
[c,m,h,nms] = multcompare(stats);

[p,t,stats] = anova1(MoveHead3DError)
[h,p,ci,stats] =ttest2(MoveHead3DError(:,1),MoveHead3DError(:,2))
[h,p,ci,stats] =ttest2(MoveHead3DError(:,3),MoveHead3DError(:,2))
[h,p,ci,stats] =ttest2(MoveHead3DError(:,1),MoveHead3DError(:,3))

[p,t,stats] = anova1(MoveReach3DError)
[h,p,ci,stats] =ttest2(MoveReach3DError(:,1),MoveReach3DError(:,2))
[h,p,ci,stats] =ttest2(MoveReach3DError(:,3),MoveReach3DError(:,2))
[h,p,ci,stats] =ttest2(MoveReach3DError(:,1),MoveReach3DError(:,3))

[p,t,stats] = anova1(PinnedReachError)
[h,p,ci,stats] =ttest2(PinnedReachError(:,1),PinnedReachError(:,2))
[h,p,ci,stats] =ttest2(PinnedReachError(:,3),PinnedReachError(:,2))
[h,p,ci,stats] =ttest2(PinnedReachError(:,1),PinnedReachError(:,3))

[h,p,ci,stats] =ttest2(Controller3DHeadError(:,1),MoveHead3DError(:,1))
[h,p,ci,stats] =ttest2(Controller3DHeadError(:,3),MoveHead3DError(:,3))
[h,p,ci,stats] =ttest2(Controller3DHeadError(:,2),MoveHead3DError(:,2))

mean(Controller3DHeadError(:,1))-mean(MoveHead3DError(:,1))
mean(Controller3DHeadError(:,3))-mean(MoveHead3DError(:,3))
mean(Controller3DHeadError(:,2))-mean(MoveHead3DError(:,2))

mean(Controller3DHeadError(:,1)-MoveHead3DError(:,1))
std(Controller3DHeadError(:,1)-MoveHead3DError(:,1))/sqrt(14)
mean(Controller3DHeadError(:,3)-MoveHead3DError(:,3))
std(Controller3DHeadError(:,3)-MoveHead3DError(:,3))/sqrt(14)

mean(MoveHead3DError(:,2)./((MoveHead3DError(:,1)+MoveHead3DError(:,3))./2))
std(MoveHead3DError(:,2)./((MoveHead3DError(:,1)+MoveHead3DError(:,3))./2))/sqrt(14)

mean(MoveHead3DError(:,2))/mean((MoveHead3DError(:,1)+MoveHead3DError(:,3))./2)

mean(MoveReach3DError(:,2)./((MoveReach3DError(:,1)+MoveReach3DError(:,3))/2))
std(MoveReach3DError(:,2)./((MoveReach3DError(:,1)+MoveReach3DError(:,3))/2))/sqrt(14)
mean(MoveReach3DError(:,2))/mean((MoveReach3DError(:,1)+MoveReach3DError(:,3))./2)

%Here we use paired ttest since we want to compare effect of condition
%within each participant
[h,p,ci,stats] =ttest((MoveHead3DError(:,2)./((MoveHead3DError(:,1)+MoveHead3DError(:,3))./2)),(MoveReach3DError(:,2)./((MoveReach3DError(:,1)+MoveReach3DError(:,3))./2)) )


% Save summary and plots
ErrorSummaryData(1,:)=PinnedReachError;
ErrorSummaryData( 1, end+1:end+3 )=MoveReach3DError;
ErrorSummaryData( 1, end+1:end+3 )=MoveHead3DError;
ErrorSummaryData( 1, end+1:end+3 )=Controller3DHeadError;
ErrorSummaryData( 1, end+1:end+3 )=PinnedReachSTE;
ErrorSummaryData( 1, end+1:end+3 )=MoveReach3DSTE;
ErrorSummearyData( 1, end+1:end+3 )=MoveHead3DSTE;
ErrorSummaryData( 1, end+1:end+3 )=Controller3DHeadSTE;

colNames = {'PinnedMeanErrorX','PinnedMeanErrorZ','PinnedMeanErrorY','MoveReachErrorX','MoveReachErrorZ','MoveReachErrorY','MoveHeadErrorX','MoveHeadErrorY','MoveHeadErrorZ','ControllerErrorX','ControllerErrorZ','ControllerErrorY','PinnedSTEX','PinnedSTEZ','PinnedSTEY','MoveReachSTEX','MoveReachSTEZ','MoveReachSTEY','MoveHeadSTEX','MoveHeadSTEY','MoveHeadSTEZ','ControllerSTEX','ControllerSTEZ','ControllerSTEY'};
rowNames = {'1'};
sTable = array2table(ErrorSummaryData,'RowNames',rowNames,'VariableNames',colNames);

name=string(s)+'\SummaryStats\ErrorSummaryData.csv';
% Save it
writetable(sTable,name);

%% Error timeplot for all participants

LfontSize= 13;
TfontSize= 14;
ticksizeX= 9;
ticksizeY= 9;
SideLineW=1.2;

CartErrorX= xlsread('CartesianErrorX.csv');
CartErrorY=xlsread('CartesianErrorY.csv');
CartErrorZ=xlsread('CartesianErrorZ.csv');

%Trim the End 
CartErrorX(:,98:end)=[];
CartErrorY(:,98:end)=[];
CartErrorZ(:,98:end)=[];

CartErrorX=abs(CartErrorX);
CartErrorY=abs(CartErrorY);
CartErrorZ=abs(CartErrorZ);
YDraw=linspace(1,length(CartErrorX),length(CartErrorX));

figure(10)
subplot(3,1,1)

shadedErrorBar(YDraw,mean(CartErrorX),std(CartErrorX))
title('Error across time for all participants','fontweight','bold','FontSize',TfontSize)
xlabel('Trial Number','fontweight','bold','FontSize',LfontSize)
ylabel({'X Absolute ','Error'},'fontweight','bold','FontSize',LfontSize)
ax = gca;
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
ax.LineWidth = SideLineW;
xlim([0,90])
ylim([-0.5,2])
subplot(3,1,2)
shadedErrorBar(YDraw,mean(CartErrorY),std(CartErrorY))

xlabel('Trial Number','fontweight','bold','FontSize',LfontSize)
ylabel({'Y Absolute ','Error'},'fontweight','bold','FontSize',LfontSize)
ax = gca;
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
ax.LineWidth = SideLineW;
xlim([0,90])
ylim([-0.5,2])
subplot(3,1,3)
shadedErrorBar(YDraw,mean(CartErrorZ),std(CartErrorZ))

xlabel('Trial Number','fontweight','bold','FontSize',LfontSize)
ylabel({'Z Absolute ','Error'},'fontweight','bold','FontSize',LfontSize)
ax = gca;
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
ax.LineWidth = SideLineW;
xlim([0,90])
ylim([-0.5,2])

s = pwd;
name=string(s)+'\TimePlotAll.jpg';
saveas(gcf,name)

s = pwd;
name=string(s)+'\TimePlotAll.svg';
saveas(gcf,name)

%% New Main plot error in dimensions with error bars

figure(11)
title('Error in dimensions for all participants')
colorBlue=[.4,.5,.8];
x0=400;
y0=200;
width=300;
height=600;
wWidth=1.8;
LfontSize= 12;
TfontSize= 16;
ticksizeX= 12;
ticksizeY= 10;
DotSize=300;
MeanSize= 25;
SideLineW=1;
Csize= 11;

set(gcf,'position',[x0,y0,width,height])

x=[.5,1.5,2.5];
location=ones(3,14);
location=x.*location'


subplot(1,4,2)
scatter(location(:,1),[MoveReach3DError(:,1)]*100,'.','SizeData',DotSize,'MarkerEdgeColor',colorBlue)   
hold on
scatter(location(:,2),[MoveReach3DError(:,3)]*100,'.','SizeData',DotSize,'MarkerEdgeColor',colorBlue)  
scatter(location(:,3),[MoveReach3DError(:,2)]*100,'.','SizeData',DotSize,'MarkerEdgeColor',colorBlue)  
xticks([.5,1.5,2.5])
xticklabels({'X','Y','Z'})
errorbar(location(1,:),[MoveReach3DMean(1),MoveReach3DMean(3),MoveReach3DMean(2)]*100, MoveReach3DSTE*100,'.','Color','k','MarkerSize',MeanSize,'LineWidth',wWidth,'CapSize',Csize); 
ylim([0,140])
xlim([0,3])
title('Reach Hybrid','fontweight','bold','FontSize',TfontSize)
xlabel('Dimension','fontweight','bold','FontSize',LfontSize)
ylabel('Error in cm','fontweight','bold','FontSize',LfontSize)
set(gca, 'YGrid', 'off', 'XGrid', 'on')
ax = gca;
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
ax.LineWidth = SideLineW;

subplot(1,4,3)
scatter(location(:,1),[MoveHead3DError(:,1)]*100,'.','SizeData',DotSize,'MarkerEdgeColor',colorBlue)   
hold on
scatter(location(:,2),[MoveHead3DError(:,3)]*100,'.','SizeData',DotSize,'MarkerEdgeColor',colorBlue)  
scatter(location(:,3),[MoveHead3DError(:,2)]*100,'.','SizeData',DotSize,'MarkerEdgeColor',colorBlue)  
xticks([.5,1.5,2.5])
xticklabels({'X','Y','Z'})
errorbar(location(1,:),[MoveHead3DMean(1),MoveHead3DMean(3),MoveHead3DMean(2)]*100, MoveHead3DSTE*100,'.','Color','k','MarkerSize',MeanSize,'LineWidth',wWidth,'CapSize',Csize); 
ylim([0,140])
xlim([0,3])
title('Hybrid','fontweight','bold','FontSize',TfontSize)
xlabel('Dimension','fontweight','bold','FontSize',LfontSize)
ylabel('Error in cm','fontweight','bold','FontSize',LfontSize)
set(gca, 'YGrid', 'off', 'XGrid', 'on')
ax = gca;
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
ax.LineWidth = SideLineW;

subplot(1,4,1)
scatter(location(:,1),[PinnedReachError(:,1)]*100,'.','SizeData',DotSize,'MarkerEdgeColor',colorBlue);  
hold on
scatter(location(:,2),[PinnedReachError(:,3)]*100,'.','SizeData',DotSize,'MarkerEdgeColor',colorBlue); 
scatter(location(:,3),[PinnedReachError(:,2)]*100,'.','SizeData',DotSize,'MarkerEdgeColor',colorBlue);
xticks([.5,1.5,2.5])
xticklabels({'X','Y','Z'})
errorbar(location(1,:),[PinnedReachMean(1),PinnedReachMean(3),PinnedReachMean(2)]*100, PinnedReachSTE*100,'.','Color','k','MarkerSize',MeanSize,'LineWidth',wWidth,'CapSize',Csize); 
ylim([0,140])
xlim([0,3])
title('Reach','fontweight','bold','FontSize',TfontSize)
xlabel('Dimension','fontweight','bold','FontSize',LfontSize)
ylabel('Error in cm','fontweight','bold','FontSize',LfontSize)
set(gca, 'YGrid', 'off', 'XGrid', 'on')
ax = gca;
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
ax.LineWidth = SideLineW;

subplot(1,4,4)
scatter(location(:,1),[Controller3DHeadError(:,1)]*100,'.','SizeData',DotSize,'MarkerEdgeColor',colorBlue)   
hold on
scatter(location(:,2),[Controller3DHeadError(:,3)]*100,'.','SizeData',DotSize,'MarkerEdgeColor',colorBlue)  
scatter(location(:,3),[Controller3DHeadError(:,2)]*100,'.','SizeData',DotSize,'MarkerEdgeColor',colorBlue)  
xticks([.5,1.5,2.5])
xticklabels({'X','Y','Z'})
errorbar(location(1,:),[Controller3DHeadMean(1),Controller3DHeadMean(3),Controller3DHeadMean(2)]*100, Controller3DHeadSTE*100,'.','Color','k','MarkerSize',MeanSize,'LineWidth',wWidth,'CapSize',Csize); 
ylim([0,140])
xlim([0,3])
title('Virtual','fontweight','bold','FontSize',TfontSize)
xlabel('Dimension','fontweight','bold','FontSize',LfontSize)
ylabel('Error in cm','fontweight','bold','FontSize',LfontSize)
set(gca, 'YGrid', 'off', 'XGrid', 'on')
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
x0=400;
y0=200;
width=800;
height=400;
set(gcf,'position',[x0,y0,width,height])
ax = gca;
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
ax.LineWidth = SideLineW;
s = pwd;
name=string(s)+'\ErrorDotPlotAll.jpg';
saveas(gcf,name)

name=string(s)+'\ErrorDotPlotAll.svg';
saveas(gcf,name)

% run pairwise ttest to see if difference is significant
[h,p,ci,stats] =ttest2(Controller3DHeadError(:,1),Controller3DHeadError(:,2))
[h,p,ci,stats] =ttest2(Controller3DHeadError(:,3),Controller3DHeadError(:,2))
% run pairwise ttest to see if x and y are the same
[h,p,ci,stats] =ttest2(Controller3DHeadError(:,1),Controller3DHeadError(:,3))


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

%% draw histograms to analyze Skewness

wWidth=2.5;
LfontSize= 11;
TfontSize= 14;
ticksizeX= 8;
ticksizeY= 8;
SideLineW=6;

CartErrorXall= xlsread('CartesianErrorX.csv');
CartErrorYall= xlsread('CartesianErrorY.csv');
CartErrorZall= xlsread('CartesianErrorZ.csv');

CartErrorXall=reshape(CartErrorXall,1,[]);
CartErrorYall=reshape(CartErrorYall,1,[]);
CartErrorZall=reshape(CartErrorZall,1,[]);

Step=.05;

figure(12)
subplot(3,1,1)
NumBins=floor((max(CartErrorXall)-min(CartErrorXall))/Step)
histfit(CartErrorXall,NumBins,'tlocationscale')
xlim([-3,3])
ylim([0,100])
title('X dimension','fontweight','bold','FontSize',TfontSize)
xlabel('Distance from veridical(m)','fontweight','bold','FontSize',LfontSize)
ylabel({'Error','frequency'},'fontweight','bold','FontSize',LfontSize)


subplot(3,1,2)
NumBins=floor((max(CartErrorYall)-min(CartErrorYall))/Step)
histfit(CartErrorYall,NumBins,'tlocationscale')
xlim([-3,3])
ylim([0,100])
title('Y dimension','fontweight','bold','FontSize',TfontSize)
xlabel('Distance from veridical(m)','fontweight','bold','FontSize',LfontSize)
ylabel({'Error','frequency'},'fontweight','bold','FontSize',LfontSize)
set(gca, 'YGrid', 'off', 'XGrid', 'on')


subplot(3,1,3)
NumBins=floor((max(CartErrorZall)-min(CartErrorZall))/Step)
histfit(CartErrorZall,NumBins,'tlocationscale')
xlim([-3,3])
ylim([0,100])
title('Z dimension','fontweight','bold','FontSize',TfontSize)
xlabel('Distance from veridical(m)','fontweight','bold','FontSize',LfontSize)
ylabel({'Error','frequency'},'fontweight','bold','FontSize',LfontSize)
set(gca, 'YGrid', 'off', 'XGrid', 'on')



s = pwd;
name=string(s)+'\SkewnessPlot.jpg';
saveas(gcf,name)
name=string(s)+'\SkewnessPlot.svg';
saveas(gcf,name)
% One measure of skewness, called Pearson's first coefficient of skewness, is to subtract the mean from the mode, and then divide this difference by the standard deviation of the data. The reason for dividing the difference is so that we have a dimensionless quantity.
% subtract the mean from the mode, and then divide this difference by the standard deviation of the data.
% When you set flag to 0, skewness corrects for the systematic bias, and the following equation applies:
PCSX=skewness(CartErrorXall,0);
PCSY=skewness(CartErrorYall,0);
PCSZ=skewness(CartErrorZall,0);

PCS(:,1)=PCSX;
PCS(:,2)=PCSY;
PCS(:,3)=PCSZ;

colNames = {'PCSX','PCSY','PCSZ'};
rowNames = {'1'};

sTable = array2table(PCS,'RowNames',rowNames,'VariableNames',colNames);
name=string(s)+'\SkewnessAll.csv';
writetable(sTable,name);

%% Skewness
SkewnessData= xlsread('AllSkewness.csv');
x0=400;
y0=200;
width=300;
height=600;
wWidth=2.3;
LfontSize= 16;
TfontSize= 17;
ticksizeX= 14;
ticksizeY= 14;
DotSize=400;
MeanSize= 28;
SideLineW=2;

figure(13)
colorBlue=[.4,.5,.8];
x0=100;
y0=100;
DotSize=170;
x=[15,25,35];
location=ones(3,14);
location=x.*location'
violins= violinplot([SkewnessData(:,1),SkewnessData(:,3),SkewnessData(:,2)]*100,'ShowMean' ,'BoxWidth' ,0.0,'Width',.3);

xticks([15,25,35])
xticklabels({'X','Y','Z'})

SkewnessSTE= std(SkewnessData)/sqrt(length(SkewnessData))

errorbar([mean(SkewnessData(:,1)),mean(SkewnessData(:,2)),mean(SkewnessData(:,3))]*100, SkewnessSTE*100,'.','Color','k','MarkerSize',MeanSize,'LineWidth',wWidth); 
ylim([-200,200])
xlim([0,4])
width=300;
height=600;
ax = gca;
set(gcf,'position',[x0,y0,width,height])
title({'Pearson Coefficient ','of Skewness '},'fontweight','bold','FontSize',TfontSize)
xlabel('Dimension','fontweight','bold','FontSize',LfontSize)
ylabel('PCS','fontweight','bold','FontSize',LfontSize)
set(gca, 'YGrid', 'off', 'XGrid', 'on')
set(get(gca, 'XAxis'), 'FontWeight', 'bold');
set(get(gca, 'YAxis'), 'FontWeight', 'bold');
ax.LineWidth = SideLineW;
grid on
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;

s = pwd;
name=string(s)+'\PCSAll1.jpg';
saveas(gcf,name)
name=string(s)+'\PCSAll1.svg';
saveas(gcf,name)

%% Then look at how the their mean affects this variance

%% movement Analysis

Move= xlsread('MoveDataAll.csv');

xmove=Move(:,1);
ymove=Move(:,3);

MoveRatio=Move(:,2)./mean([xmove,ymove],2);
[h,p]=ttest(Move(:,8),mean([xmove,ymove],2))

Fit= xlsread('SR.csv');

Verts= (Fit(:,3)+Fit(:,4))./2;
Horts= (Fit(:,2));

probRatio= (Verts./Horts)

DotSize=600;
colorBlue=[.4,.4,.8];

figure(19)
scatter(MoveRatio,probRatio,'.','SizeData',DotSize,'MarkerEdgeColor',colorBlue);  

grid on
xlabel("Ratio of Normalized Vertical to Horizontal Movement")
ylabel("Ratio of Vertical to Horizontal Transition Probabilities")

grid on
hold on
p = polyfit( MoveRatio,probRatio,1);
x1 = linspace(0,5,2);
y1 = polyval(p,x1);

plot(x1,y1,'LineWidth',.5,'color','r')

ylim([0.5 2])
xlim([0.75 2.0])

[R,P] = corrcoef(MoveRatio,probRatio)
title("Correlation movement stats with SR probabilities ")

str=['  r = ',num2str(R(1,2))]
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

s = pwd;
name=string(s)+'\CorrMoveSRNewVersion.jpg';
saveas(gcf,name)
name=string(s)+'\CorrMoveSRNewVersion.svg';
saveas(gcf,name)


%% Gamma Analysis

Move= xlsread('MoveDataAll.csv');
Allmove=Move(:,1)+Move(:,2)+Move(:,3);
Gammas= Fit(:,1);
DotSize=600;
colorBlue=[.4,.4,.8];

figure(25)
scatter((Allmove),Gammas,'.','SizeData',DotSize,'MarkerEdgeColor',colorBlue);
hold on
[fitresult, gof] = FitPower2(Allmove, Gammas);
sqrt(gof.rsquare)
help FitPower2

str=['  R2 ',num2str(gof.rsquare)]
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

xlabel("Sum of Movement")
ylabel("Gamma from Fit")

grid on
hold on

ylim([0.89 1.01])
xlim([200 1400])

s = pwd;
name=string(s)+'\GammaCorrelation.jpg';
saveas(gcf,name)
name=string(s)+'\GammaCorrelation.svg';
saveas(gcf,name)

ErrorData= xlsread('SummaryErrors.csv');
allError=sum(ErrorData(:,11:13),2);
figure(26)
scatter(allError,Gammas,'.','SizeData',DotSize,'MarkerEdgeColor',colorBlue);


%% Movemelt Plot

Move= xlsread('MoveDataAll.csv');
MoveMean=mean(Move);

MoveMean(:,[2 3]) = MoveMean(:,[3 2]);
MoveMean(:,[5 6]) = MoveMean(:,[6 5]);
MoveMean(:,[8 9]) = MoveMean(:,[9 8]);
MoveMean(:,[11 12]) = MoveMean(:,[12 11]);
MoveMean(:,[14 15]) = MoveMean(:,[15 14]);

toPlot(1,:)=MoveMean(1:3);
toPlot(2,:)=MoveMean(4:6);
toPlot(3,:)=MoveMean(7:9); 
toPlot(4,:)=MoveMean(10:12); 
toPlot(5,:)=MoveMean(13:15); 
toPlot(6,:)=MoveMean(16:18);
MoveSTE(1,1)=(std(Move(:,1))/(sqrt(length(Move(:,1)))));
MoveSTE(1,2)=(std(Move(:,3))/(sqrt(length(Move(:,3)))));
MoveSTE(1,3)=(std(Move(:,2))/(sqrt(length(Move(:,2)))));

MoveSTE(2,1)=(std(Move(:,4))/(sqrt(length(Move(:,1)))));
MoveSTE(2,2)=(std(Move(:,6))/(sqrt(length(Move(:,1)))));
MoveSTE(2,3)=(std(Move(:,5))/(sqrt(length(Move(:,1)))));

MoveSTE(3,1)=(std(Move(:,7))/(sqrt(length(Move(:,1)))));
MoveSTE(3,2)=(std(Move(:,9))/(sqrt(length(Move(:,1)))));
MoveSTE(3,3)=(std(Move(:,8))/(sqrt(length(Move(:,1)))));

MoveSTE(4,1)=(std(Move(:,10))/(sqrt(length(Move(:,1)))));
MoveSTE(4,2)=(std(Move(:,12))/(sqrt(length(Move(:,1)))));
MoveSTE(4,3)=(std(Move(:,11))/(sqrt(length(Move(:,1)))));

MoveSTE(5,1)=(std(Move(:,12))/(sqrt(length(Move(:,1)))));
MoveSTE(5,2)=(std(Move(:,14))/(sqrt(length(Move(:,1)))));
MoveSTE(5,3)=(std(Move(:,13))/(sqrt(length(Move(:,1)))));

MoveSTE(6,1)=(std(Move(:,14))/(sqrt(length(Move(:,1)))));
MoveSTE(6,2)=(std(Move(:,16))/(sqrt(length(Move(:,1)))));
MoveSTE(6,3)=(std(Move(:,15))/(sqrt(length(Move(:,1)))));

figure(22)
subplot(1,5,1)
violins= violinplot([Move(:,1),Move(:,3),Move(:,2)],'ShowMean' ,'MedianColor',[0,0,0],'BoxWidth' ,0.0,'Width',.2);
hold on
errorbar(toPlot(1,:), MoveSTE(1,:),'.','Color','k','MarkerSize',9,'linewidth',.5); 
grid on
xticklabels({'X','Y','Z'})
ylim([-100,500])
xlim([0,4])
yticks(-100:100:500);
title('All Movement')

subplot(1,5,2)
violins= violinplot([Move(:,4),Move(:,6),Move(:,5)],'ShowMean' ,'MedianColor',[0,0,0],'BoxWidth' ,0.0,'Width',.2);
hold on
errorbar(toPlot(2,:), MoveSTE(2,:),'.','Color','k','MarkerSize',9,'linewidth',.5); 
grid on
xticklabels({'X','Y','Z'})
ylim([-100,500])
xlim([0,4])
yticks(-100:100:500);
title('Search')

subplot(1,5,3)
violins= violinplot([Move(:,7),Move(:,9),Move(:,8)],'ShowMean' ,'MedianColor',[0,0,0],'BoxWidth' ,0.0,'Width',.2);
hold on
errorbar(toPlot(3,:), MoveSTE(3,:),'.','Color','k','MarkerSize',9,'linewidth',.5); 
grid on
xticklabels({'X','Y','Z'})
ylim([-100,500])
xlim([0,4])
yticks(-100:100:500);
title('Navigate')

subplot(1,5,4)
violins= violinplot([Move(:,10),Move(:,12),Move(:,11)],'ShowMean' ,'MedianColor',[0,0,0],'BoxWidth' ,0.0,'Width',.2);
hold on
errorbar(toPlot(4,:), MoveSTE(4,:),'.','Color','k','MarkerSize',9,'linewidth',.5); 
grid on
xticklabels({'X','Y','Z'})
ylim([-75,75])
xlim([0,4])
yticks(-100:100:500);
title('Search Abs')

subplot(1,5,5)
violins= violinplot([Move(:,13),Move(:,15),Move(:,14)],'ShowMean' ,'MedianColor',[0,0,0],'BoxWidth' ,0.0,'Width',.2);
hold on
errorbar(toPlot(5,:), MoveSTE(5,:),'.','Color','k','MarkerSize',9,'linewidth',.5); 
grid on
xticklabels({'X','Y','Z'})
ylim([-75,75])
xlim([0,4])
yticks(-100:100:500);
title('Navigate Abs')

set(gcf,'Position',[100 100 800 400])
s = pwd;
name=string(s)+'\MovementAnalysisNew1.jpg';
saveas(gcf,name)
name=string(s)+'\MovementAnalysisNew1.svg';
saveas(gcf,name)

