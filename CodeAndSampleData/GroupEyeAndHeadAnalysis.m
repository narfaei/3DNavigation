%% Group Eye and Head Analysis

%%  Housekeeping
clear
clc
fclose('all');
close all

%% Color Pallettes

myBlue=[0,0.447000000000000,0.741000000000000];
myCyan=[0.301000000000000,0.745000000000000,0.933000000000000];
myPurple=[0.494000000000000,0.184000000000000,0.556000000000000];
myOrange=[0.635000000000000,0.0980000000000000,0.184000000000000];
myRed=[0.635000000000000,0.0780000000000000,0.184000000000000];
myGreen= [0.466000000000000,0.674000000000000,0.188000000000000];

mycolormap = customcolormap([0 0.7 0.9 1], [myRed; myPurple;myCyan; myBlue ]);
set(0, 'defaultFigureRenderer', 'painters')

wWidth=2.5;
LfontSize= 22;
TfontSize= 24;
ticksizeX= 14;
ticksizeY= 14;
ticksizeZ= 14;
DotSize=20;

%% Load data
CovGazePosNav= xlsread('CovGazePosNav.csv');
CovGazePosSearch= xlsread('CovGazePosSearch.csv');

%% Plot the group gaze position cov matrices
CovGazePosNavMean=mean(CovGazePosNav);
reorganized=CovGazePosNavMean(1:3);
reorganized(2,1:3)=CovGazePosNavMean(4:6);
reorganized(3,1:3)=CovGazePosNavMean(7:9);

figure()
x0=400;
y0=200;
width=700;
height=220;
set(gcf,'position',[x0,y0,width,height])
subplot(1,2,1)
imagesc(reorganized(1:2,1:2));

set(gca, 'XTick', 1:3); % center x-axis ticks on bins
set(gca, 'YTick', 1:3); % center y-axis ticks on bins
set(gca, 'XTickLabel', {'X' ,'Y' ,'Z'}); % set x-axis labels
set(gca, 'YTickLabel',  {'X' ,'Y' ,'Z'}); % set y-axis labels
title('Covariance of GazePosNav', 'FontSize', 10); % set title
caxis manual
caxis([0 0.05]);
colorbar;
colormap(mycolormap);% Choose jet or any other color scheme
colorbar;

CovGazePosSearchMean=mean(CovGazePosSearch);
reorganized=CovGazePosSearchMean(1:3);
reorganized(2,1:3)=CovGazePosSearchMean(4:6);
reorganized(3,1:3)=CovGazePosSearchMean(7:9);

subplot(1,2,2)
imagesc(reorganized(1:2,1:2));
set(gca, 'XTick', 1:3); % center x-axis ticks on bins
set(gca, 'YTick', 1:3); % center y-axis ticks on bins
set(gca, 'XTickLabel', {'X' ,'Y' ,'Z'}); % set x-axis labels
set(gca, 'YTickLabel',  {'X' ,'Y' ,'Z'}); % set y-axis labels
title('Covariance of GazePosSearch', 'FontSize', 10); % set title
caxis manual
caxis([0 0.05]);
colormap(mycolormap);% Choose jet or any other color scheme
colorbar;

s = pwd;
name=string(s)+'\CovGazePos.jpg';
saveas(gcf,name)

%% Plot the group gaze movement cov matrices

CovGazeMoveNav= xlsread('CovGazeMoveNav.csv');
CovGazeMoveSearch= xlsread('CovGazeMoveSearch.csv');

CovGazeMoveNavMean=mean(CovGazeMoveNav);
reorganized=CovGazeMoveNavMean(1:3);
reorganized(2,1:3)=CovGazeMoveNavMean(4:6);
reorganized(3,1:3)=CovGazeMoveNavMean(7:9);

figure()
x0=400;
y0=200;
width=700;
height=220;
set(gcf,'position',[x0,y0,width,height])
subplot(1,2,1)
imagesc(reorganized(1:2,1:2));

set(gca, 'XTick', 1:3); % center x-axis ticks on bins
set(gca, 'YTick', 1:3); % center y-axis ticks on bins
set(gca, 'XTickLabel', {'X' ,'Y' ,'Z'}); % set x-axis labels
set(gca, 'YTickLabel',  {'X' ,'Y' ,'Z'}); % set y-axis labels
title('Covariance of GazeMoveNav', 'FontSize', 10); % set title
caxis manual
caxis([0 0.0005]);
colormap(mycolormap);% Choose jet or any other color scheme
colorbar;

CovGazeMoveSearchMean=mean(CovGazeMoveSearch);
reorganized=CovGazeMoveSearchMean(1:3)
reorganized(2,1:3)=CovGazeMoveSearchMean(4:6)
reorganized(3,1:3)=CovGazeMoveSearchMean(7:9)

subplot(1,2,2)
imagesc(reorganized(1:2,1:2));

set(gca, 'XTick', 1:3); % center x-axis ticks on bins
set(gca, 'YTick', 1:3); % center y-axis ticks on bins
set(gca, 'XTickLabel', {'X' ,'Y' ,'Z'}); % set x-axis labels
set(gca, 'YTickLabel',  {'X' ,'Y' ,'Z'}); % set y-axis labels
title('Covariance of GazeMoveSearch', 'FontSize', 10); % set title
caxis manual
caxis([0 0.0005]);
colormap(mycolormap);% Choose jet or any other color scheme

colorbar;

s = pwd;
name=string(s)+'\CovGazeMove.jpg';
saveas(gcf,name)
%% Cov of gaze movement during epochs

headPosFactors= xlsread('headPosFactors.csv');

headPosFactorsMean=nanmean(headPosFactors);
reorganized=CovGazeMoveNavMean(1:3)
reorganized(2,1:3)=CovGazeMoveNavMean(4:6)
reorganized(3,1:3)=CovGazeMoveNavMean(7:9)

figure()
x0=400;
y0=200;
width=700;
height=220;
set(gcf,'position',[x0,y0,width,height])
subplot(1,2,1)
imagesc(reorganized(1:2,1:2));

set(gca, 'XTick', 1:3); % center x-axis ticks on bins
set(gca, 'YTick', 1:3); % center y-axis ticks on bins
set(gca, 'XTickLabel', {'X' ,'Y' ,'Z'}); % set x-axis labels
set(gca, 'YTickLabel',  {'X' ,'Y' ,'Z'}); % set y-axis labels
title('Covariance of GazeMoveNav', 'FontSize', 10); % set title
caxis manual
caxis([0 0.0005]);
colormap(mycolormap);% Choose jet or any other color scheme
colorbar;

CovGazeMoveSearch= xlsread('CovGazeMoveSearch.csv');

CovGazeMoveSearchMean=mean(CovGazeMoveSearch);
reorganized=CovGazeMoveSearchMean(1:3)
reorganized(2,1:3)=CovGazeMoveSearchMean(4:6)
reorganized(3,1:3)=CovGazeMoveSearchMean(7:9)


subplot(1,2,2)
imagesc(reorganized(1:2,1:2));

set(gca, 'XTick', 1:3); % center x-axis ticks on bins
set(gca, 'YTick', 1:3); % center y-axis ticks on bins
set(gca, 'XTickLabel', {'X' ,'Y' ,'Z'}); % set x-axis labels
set(gca, 'YTickLabel',  {'X' ,'Y' ,'Z'}); % set y-axis labels
title('Covariance of GazeMoveSearch', 'FontSize', 10); % set title
caxis manual
caxis([0 0.0005]);
colormap(mycolormap);% Choose jet or any other color scheme
colormap(mycolormap);% Choose jet or any other color scheme
colorbar;

s = pwd;
name=string(s)+'\CovGazeMove.jpg';
saveas(gcf,name)

%% Plot Group Heatmaps

NavEyeData= xlsread('NavEyePos.csv');
mycolormap = customcolormap([0 0.8 0.95 1], [myRed; myPurple;myCyan; myBlue]);
pts = linspace(-1, 1, 41);

figure()
x0=50; y0=50; width=800; height=700; set(gcf,'position',[x0,y0,width,height]);

subplot(2,2,2)
N = histcounts2(NavEyeData(:,2), NavEyeData(:,1), pts, pts);
N=N./sum(sum(N));
DensitiesNavEyePos=N;
imagesc(pts, pts, N);
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
% scatter(GazePos(1,CurrEvent=="inTrial"),GazePos(2,CurrEvent=="inTrial"),'.','SizeData',DotSize,'MarkerEdgeColor','b')
shading interp
axis equal
xlabel('X','fontweight','bold','FontSize',LfontSize)
ylabel('Y','fontweight','bold','FontSize',LfontSize)
title(' During Navigation','fontweight','bold','FontSize',TfontSize)

ylim([-1,1])
xlim([-1,1])
xticks(-1:0.25:1);
yticks(-1:0.25:1);
colormap(mycolormap);
colorbar;
caxis([0 0.03])

subplot(2,2,1)
searchEyePos= xlsread('SearchEyePos.csv');
pts = linspace(-1, 1, 41);
N = histcounts2(searchEyePos(:,2), searchEyePos(:,1), pts, pts);
N=N./sum(sum(N));
DensitiesSearchEyePos=N;
imagesc(pts, pts, N);
axis equal;
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
axis equal;
xlabel('X','fontweight','bold','FontSize',LfontSize)
ylabel('Y','fontweight','bold','FontSize',LfontSize)
title('During Search','fontweight','bold','FontSize',TfontSize)
ylim([-1,1]);xlim([-1,1]);
xticks(-1:0.25:1); yticks(-1:0.25:1);
colorbar;
caxis([0 0.03])

subplot(2,2,3)
ReachsearchEyePos= xlsread('ReachSearchEyePosition.csv');
pts = linspace(-1, 1, 41);
N = histcounts2(ReachsearchEyePos(:,2), ReachsearchEyePos(:,1), pts, pts);
N=N./sum(sum(N));
DensitiesReachsearchEyePos=N;
imagesc(pts, pts, N);
axis equal;
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
axis equal;
xlabel('X','fontweight','bold','FontSize',LfontSize)
ylabel('Y','fontweight','bold','FontSize',LfontSize)
title('During Search-Reach','fontweight','bold','FontSize',TfontSize)
ylim([-1,1]);xlim([-1,1]);
xticks(-1:0.25:1); yticks(-1:0.25:1);
colorbar;
caxis([0 0.03])

subplot(2,2,4)
ReachNavEyeData= xlsread('ReachEyePositionNav.csv');
pts = linspace(-1, 1, 41);
N = histcounts2(ReachNavEyeData(:,2), ReachNavEyeData(:,1), pts, pts);
N=N./sum(sum(N));
DensitiesReachNavEyeData=N;
imagesc(pts, pts, N);
axis equal;
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
axis equal;
xlabel('X','fontweight','bold','FontSize',LfontSize)
ylabel('Y','fontweight','bold','FontSize',LfontSize)
title('During Locate-Reach','fontweight','bold','FontSize',TfontSize)
ylim([-1,1]);xlim([-1,1]);
xticks(-1:0.25:1); yticks(-1:0.25:1);
colorbar;
caxis([0 0.03])

% Save it
s= pwd;
name=string(s)+'\BothGazePosHeatmapFinal.jpg';
saveas(gcf,name)
name=string(s)+'\BothGazePosHeatmapFinal.svg';
saveas(gcf,name)

%% Gaze Location Densities plot
pts = linspace(-1, 1, 90);
figure()
subplot(2,2,1)
plot(sum(DensitiesSearchEyePos,2))
hold on
plot(sum(DensitiesSearchEyePos,1)')
ylim([0,0.4])
title('Search-Eyemove')

subplot(2,2,2)
plot(sum(DensitiesNavEyePos,2))
hold on
plot(sum(DensitiesNavEyePos,1)')
title('Nav-Eyemove')
ylim([0,0.4])

subplot(2,2,3)
plot(sum(DensitiesReachsearchEyePos,2))
hold on
plot(sum(DensitiesReachsearchEyePos,1)')
title('Reach-Search-Eyemove')
ylim([0,0.4])

subplot(2,2,4)
plot(sum(DensitiesReachNavEyeData,2))
hold on
plot(sum(DensitiesReachNavEyeData,1)')
title('Reach-Nav-Eyemove')
ylim([0,0.4])

% Save it
s = pwd;
name=string(s)+'\GazeLocationDensitiesAll.jpg';
saveas(gcf,name)

name=string(s)+'\GazeLocationDensitiesAll.svg';
saveas(gcf,name)

%% Plot Heatmaps of eye movement

TfontSize=18;
NavEyeMove= xlsread('NavEyeMove.csv');

nanmean(abs(NavEyeMove))
NavEyeMove= movmean(NavEyeMove,10).*90;
NavEyeMove = NavEyeMove(any(NavEyeMove,2),:);
mycolormap = customcolormap([0 0.2 0.83 0.95 1], [myOrange; myRed; myPurple;myCyan; myBlue; ]);

pts = linspace(-5,5, 501); N = histcounts2(NavEyeMove(:,1), NavEyeMove(:,2), pts, pts);
T= (N); T(T== -Inf) = 0; T(T== 0) = 0; T=T./sum(sum(T));

figure()
x0=50; y0=50; width=800; height=700; set(gcf,'position',[x0,y0,width,height]);

subplot(2,2,2)
imagesc(pts, pts, T);
colormap(mycolormap); colorbar;
caxis([0 0.005]);
xlabel('X','fontweight','bold','FontSize',LfontSize);
ylabel('Y','fontweight','bold','FontSize',LfontSize);
title('Navigation-Eyemove','fontweight','bold','FontSize',TfontSize);
axis equal;
ylim([-1,1]); xlim([-1,1]);
xticks(-1:0.25:1);yticks(-1:0.25:1);

searchEyemove= xlsread('searchEyemove.csv');
nanmean(abs(searchEyemove))
searchEyemove= movmean(searchEyemove,10).*90;
searchEyemove = searchEyemove(any(searchEyemove,2),:);
N = histcounts2(searchEyemove(:,2), searchEyemove(:,1), pts, pts);
T= (N); T(T== -Inf) = 0;T=T./sum(sum(T));

subplot(2,2,1)
imagesc( pts, pts, T);
colormap(mycolormap); colorbar;
caxis([0 0.005]);
LfontSize= 12;

xlabel('X','fontweight','bold','FontSize',LfontSize)
ylabel('Y','fontweight','bold','FontSize',LfontSize)
title('Search-Eyemove','fontweight','bold','FontSize',TfontSize)
axis equal;
ylim([-1,1]); xlim([-1,1]);
xticks(-1:0.25:1);yticks(-1:0.25:1);

%
ReachNavEyeMovement= xlsread('ReachNavEyeMovement.csv');
nanmean(abs(ReachNavEyeMovement))
ReachNavEyeMovement= movmean(ReachNavEyeMovement,10).*90;
ReachNavEyeMovement = ReachNavEyeMovement(any(ReachNavEyeMovement,2),:);
N = histcounts2(ReachNavEyeMovement(:,2), ReachNavEyeMovement(:,1), pts, pts);
T= (N); T(T== -Inf) = 0;T=T./sum(sum(T));

subplot(2,2,4)
imagesc(pts, pts, T);
colormap(mycolormap); colorbar;
caxis([0 0.005]);
xlabel('X','fontweight','bold','FontSize',LfontSize)
ylabel('Y','fontweight','bold','FontSize',LfontSize)
title('Navigation-Eyemove-reach','fontweight','bold','FontSize',TfontSize)
axis equal;

ylim([-1,1]); xlim([-1,1]);
xticks(-1:0.25:1);yticks(-1:0.25:1);

%
ReachsearchEyemovement= xlsread('ReachsearchEyemovement.csv');
nanmean(abs(ReachsearchEyemovement))
ReachsearchEyemovement= movmean(ReachsearchEyemovement,10).*90;
ReachsearchEyemovement = ReachsearchEyemovement(any(ReachsearchEyemovement,2),:);
N = histcounts2(ReachsearchEyemovement(:,2), ReachsearchEyemovement(:,1), pts, pts);
T= (N)
T(T== -Inf) = 0;
T=T./sum(sum(T));

subplot(2,2,3)
imagesc( pts, pts, T);
colormap(mycolormap); colorbar;
caxis([0 0.005]);
LfontSize= 12;

xlabel('X','fontweight','bold','FontSize',LfontSize)
ylabel('Y','fontweight','bold','FontSize',LfontSize)
title('Search-Eyemove-Reach','fontweight','bold','FontSize',TfontSize)
axis equal;
ylim([-1,1]); xlim([-1,1]);
xticks(-1:0.25:1);yticks(-1:0.25:1);
% Save it
s = pwd;
name=string(s)+'\GazeMoveHeatmapAll.jpg';
saveas(gcf,name)

name=string(s)+'\GazeMoveHeatmapAll.svg';
saveas(gcf,name)

%% Timecourse of all gaze movements

figure()
subplot(4,1,1)
plot(searchEyemove(:,1))
hold on
plot(searchEyemove(:,2))
subplot(4,1,2)
plot(NavEyeMove(:,1))
hold on
plot(NavEyeMove(:,2))
subplot(4,1,3)
plot(ReachsearchEyemovement(:,1))
hold on
plot(ReachsearchEyemovement(:,2))
subplot(4,1,4)
plot(ReachNavEyeMovement(:,1))
hold on
plot(ReachNavEyeMovement(:,2))
s = pwd;
name=string(s)+'\GazeMoveRaw.jpg';
saveas(gcf,name)

name=string(s)+'\GazeMoveRaw.svg';
saveas(gcf,name)

%% Density of gaze movements
pts = linspace(-5, 5, 101);

figure()
subplot(2,2,1)
ksdensity(searchEyemove(:,2),pts)
ylim([0,5])
xlim([-2,2])
title('Search-Eyemove')
hold on
ksdensity(searchEyemove(:,1),pts)

subplot(2,2,2)
ksdensity(searchEyemove(:,2),pts)
title('Nav-Eyemove')
hold on
ksdensity(NavEyeMove(:,1),pts)
ylim([0,5])
xlim([-2,2])
subplot(2,2,3)
ksdensity(ReachsearchEyemovement(:,2),pts)
title('Reach-Search-Eyemove')
hold on
ksdensity(ReachsearchEyemovement(:,1),pts)
ylim([0,5])
xlim([-2,2])

subplot(2,2,4)
ksdensity(ReachNavEyeMovement(:,2),pts)
title('Reach-Nav-Eyemove')
hold on
ksdensity(ReachNavEyeMovement(:,1),pts)
ylim([0,5])
xlim([-2,2])
s = pwd;
name=string(s)+'\GazeMoveDensitiesAll.jpg';
saveas(gcf,name)

name=string(s)+'\GazeMoveDensitiesAll.svg';
saveas(gcf,name)

%%  Head rotation position heatmap

CovGazePosNav= xlsread('ReachheadPosNav.csv');

NavRotDiff=diff(CovGazePosNav);
NavRotDiff=rmoutliers( NavRotDiff);
meanHeadRotNav=sum(sqrt(sum(abs(NavRotDiff.^2),2)))/(100*pi*2);
rad2deg(meanHeadRotNav);
CovGazePosSearch =xlsread('ReachheadPosSearch.csv');
SearchRotDiff=diff(CovGazePosSearch);
SearchRotDiff=rmoutliers( SearchRotDiff);
meanHeadRotSearch=sum(sqrt(sum(abs(SearchRotDiff.^2),2)))/(100*pi*2)
rad2deg(meanHeadRotSearch)

[x,y,z] = sphere;
pts = linspace(-1, 1, 22);
[azimuth,elevation,r] = cart2sph(CovGazePosSearch(:,1),CovGazePosSearch(:,3),CovGazePosSearch(:,2));
N = histcounts2(azimuth,elevation, pts, pts);
N(11,16)=1400;

figure()
subplot(1,2,1)
imagesc(nanmean(N,2)./sum(nanmean(N,2)))
colorbar;
colormap(mycolormap);
caxis([0.01 0.08])

subplot(1,2,2)
imagesc(mean(N,1)./sum(mean(N,1)))
colorbar;

colormap(mycolormap);
caxis([0.01 0.08]);


%%
CovGazePosSearch= xlsread('ReachheadPosSearch.csv');
[x,y,z] = sphere;
pts = linspace(-360, 360, 100);
[azimuth,elevation,r] = cart2sph(CovGazePosSearch(:,1),CovGazePosSearch(:,3),CovGazePosSearch(:,2));
normFactor=sum(abs(diff(azimuth)))+sum(abs(diff(elevation)));
figure()
subplot(2,2,1)
ksdensity(rad2deg(diff(azimuth)).*90,pts)
ylim([0,0.04]); xlim([-360,360]);
xlabel('Azimuth-ReachSearch','fontweight','bold','FontSize',LfontSize)
subplot(2,2,3)
ksdensity(rad2deg(diff(elevation)).*90,pts)
ylim([0,0.04]); xlim([-360,360]);
xlabel('Elevation-ReachSearch','fontweight','bold','FontSize',LfontSize)

CovGazePosNav= xlsread('ReachheadPosNav.csv');
[azimuth,elevation,r] = cart2sph(CovGazePosNav(:,1),CovGazePosNav(:,3),CovGazePosNav(:,2));

subplot(2,2,2)
ksdensity(rad2deg(diff(azimuth)).*90,pts)
ylim([0,0.04]); xlim([-360,360]);
xlabel('Azimuth-ReachNav','fontweight','bold','FontSize',LfontSize)

subplot(2,2,4)
ksdensity(rad2deg(diff(elevation)).*90,pts)
ylim([0,0.04]); xlim([-360,360]);
xlabel('Elevation-ReachNav','fontweight','bold','FontSize',LfontSize)

s=pwd;
name=string(s)+'\ReachHeadDensities.jpg';
saveas(gcf,name)
name=string(s)+'\ReachHeadDensities.svg';
saveas(gcf,name)

N = histcounts2(azimuth,elevation, pts, pts);
N(11,16)=1400


%% Eye Velocity of all

myVelocity= xlsread('EyeVelocityEpoch.csv');
times= linspace(-5,10,151);
YDraw=linspace(1,length(myVelocity),length(myVelocity));

figure()
shadedErrorBar(times,mean(myVelocity)*10*9,std(myVelocity)*10*60/sqrt(size(myVelocity,1)))
hold on
plot(times,mean(myVelocity)*10*9,'color',myBlue, 'linewidth',2)
this=mean(myVelocity(:,times>0 & times<9.6))
that= mean(myVelocity(:,times<0 ))
plot(times(times>0 & times<9.6),this*10*9,'color',myCyan, 'linewidth',2)

xlim([-5,8]);
xline(0,'--','color',myRed, 'linewidth',2)
title('Eye Velocity during epochs','FontSize',TfontSize)
xlabel('Time(s)','fontweight','bold','FontSize',LfontSize)
ylabel('Average Movement(m/s)','fontweight','bold','FontSize',LfontSize)
width=500;
height=400;
x0=100;
y0=100;
set(gcf,'position',[x0,y0,width,height])

xticklabels({'Start','Detect','Locate'})
ax = gca;
s=pwd;
ax.XTick = [-5 0 7];
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;

name=string(s)+'\EyeVelocity.jpg';
saveas(gcf,name)
name=string(s)+'\EyeVelocity.svg';
saveas(gcf,name)

% Its stats
thisNav=mean(myVelocity(:,times>0 & times<5));
thisSearch= mean(myVelocity(:,times<0 & times>-5 ));

nanmean(thisNav)*10*9;
nanstd(thisNav*10*9);

nanmean(thisSearch)*10*9;
nanstd(thisSearch*10*9);

[h,p]=ttest2(thisNav*10*9,thisSearch*10*9);

RightAfter=mean(myVelocity(:,times>0.2 & times<0.5));
RightBefore= mean(myVelocity(:,times<-0.2 & times>-0.5 ));

mean(RightAfter*10*9);
mean(RightBefore*10*9);

%% Pupil analysis

PupilLeft= xlsread('PupilLeft.csv');
PupilRight=xlsread('PupilRight.csv');
nanmean(PupilLeft,2);

PupilLeft(PupilLeft<1)=NaN;
PupilRight(PupilRight<1)=NaN;

figure()
shadedErrorBar(times,nanmean(PupilLeft,2),nanstd(PupilLeft')./sqrt(14))
hold on
plot(times,nanmean(PupilLeft,2),'color',myBlue, 'linewidth',2)

shadedErrorBar(times,nanmean(PupilRight,2),nanstd(PupilRight')./sqrt(14))
hold on
plot(times,nanmean(PupilRight,2),'color',myCyan, 'linewidth',2)

xlim([-5,8]);
ylim([2.8,3.3]);
xline(0,'--','color',myRed, 'linewidth',2)
title('Pupil dilation during epochs')
xlabel('Time(s)','fontweight','bold','FontSize',LfontSize)
ylabel('Diameter (mm)','fontweight','bold','FontSize',LfontSize)
width=500;
height=400;
x0=100;
y0=100;
set(gcf,'position',[x0,y0,width,height])
xticklabels({'Start','Detect','Locate'})
ax = gca;
s=pwd;
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
ax.XTick = [-5 0 7];

name=string(s)+'\PupilDilation.jpg';
saveas(gcf,name)
name=string(s)+'\PupilDilation.svg';
saveas(gcf,name)

% Its stats
LeftNav=(PupilLeft(times>0 & times<5,:));
LeftSearch= (PupilLeft(times<0 & times>-5,: ));

[h,p]=ttest(mean(LeftNav,2),mean(LeftSearch,2));

RightAfter=mean(myVelocity(:,times>0.2 & times<0.5));
RightBefore= mean(myVelocity(:,times<-0.2 & times>-0.5 ));

mean(RightAfter*10*9)
mean(RightBefore*10*9)

