clear
clc
fclose('all');
close all

%% Load Data
Controller3DHeadContin= readtable('ContinousData-4-Controller3DHead.csv','Format','auto');

%% Define Colors

myBlue=[0,0.447000000000000,0.741000000000000]
myCyan=[0.301000000000000,0.745000000000000,0.933000000000000]
myPurple=[0.494000000000000,0.184000000000000,0.556000000000000]
myOrange=[0.850000000000000,0.325000000000000,0.0980000000000000]
myRed=[0.635000000000000,0.0780000000000000,0.184000000000000]
myGreen= [0.466000000000000,0.674000000000000,0.188000000000000]
set(0, 'defaultFigureRenderer', 'painters')

wWidth=2.5;
LfontSize= 22;
TfontSize= 24;
ticksizeX= 14;
ticksizeY= 14;
ticksizeZ= 14;
DotSize=20;
%% Extract EyeStats
Controller3DHeadContin=table2array(Controller3DHeadContin);
trialHitpoints=Controller3DHeadContin(20:22,:);
CurrEvent=Controller3DHeadContin(5,:);
eyeOpenness=Controller3DHeadContin(38:39,:);

%%  Headpositions During different Epochs

GazePos=Controller3DHeadContin(12:14,:);
GazePos=cellfun(@str2double,GazePos)

rot=Controller3DHeadContin(9:11,:);
rot=cellfun(@str2double,rot);
[rotazimuth,rotelevation,r] = cart2sph(rot(1,:),rot(3,:),rot(2,:));

trialNum=Controller3DHeadContin(4,:);
trialNum=cellfun(@str2double,trialNum);

horzFactorSearch= (diff(rotazimuth(CurrEvent=="StimShowing")));
vertFactorSearch= (diff(rotelevation(CurrEvent=="StimShowing")));

horzFactorNav= (diff(rotazimuth(CurrEvent=="inTrial")));
vertFactorNav= (diff(rotelevation(CurrEvent=="inTrial")));
s=pwd;
headposFactors= [nanstd(horzFactorNav),nanstd(vertFactorNav),nanstd(horzFactorSearch),nanstd(vertFactorSearch)];
colNames = {'horzFactorNav','vertFactorNav','horzFactorSearch','vertFactorSearch'};
sTable = array2table(headposFactors,'VariableNames',colNames);
name=string(s)+'\SummaryStats\headposFactors.csv';
writetable(sTable,name);

headposSearch= [horzFactorSearch;vertFactorSearch;rot(2,CurrEvent=="StimShowing")];
sTable = array2table(headposSearch);
name=string(s)+'\SummaryStats\headposSearch.csv';
writetable(sTable,name);

headposNav= [rot(1,CurrEvent=="inTrial");rot(3,CurrEvent=="inTrial");rot(2,CurrEvent=="inTrial")];
sTable = array2table(headposNav);
name=string(s)+'\SummaryStats\headposNav.csv';
writetable(sTable,name);
%% Eye Movements During Epochs

GazePos(GazePos==-1)=NaN;
GazeDiff= diff(GazePos');
GazeDiff=GazeDiff';

s=pwd;
navEyemove=(GazeDiff(:,CurrEvent=="inTrial")');
nanvar(navEyemove)
colNames = {'X','Y','Z'};
sTable = array2table(navEyemove,'VariableNames',colNames);
name=string(s)+'\SummaryStats\ReahcNavEyeMove.csv';
writetable(sTable,name);
searchEyemove=(GazeDiff(:,CurrEvent=="StimShowing")');
nanvar(searchEyemove)
colNames = {'X','Y','Z'};
sTable = array2table(searchEyemove,'ReachVariableNames',colNames);
name=string(s)+'\SummaryStats\searchEyemove.csv';
writetable(sTable,name);

%%

Controller3DHeadContin=table2array(Controller3DHeadContin);
trialHitpoints=Controller3DHeadContin(20:22,:);
CurrEvent=Controller3DHeadContin(5,:);

eyeOpenness=Controller3DHeadContin(38:39,:);
trialTime=cellfun(@str2double,Controller3DHeadContin(3,:));

rot=Controller3DHeadContin(9:11,:);
rot=cellfun(@str2double,rot);

[azimuth,elevation,r] =  cart2sph(rot(1,:),rot(3,:),rot(2,:));

nanmean(azimuth)*90;
nanmean(elevation)*90;

AzDiff= diff(azimuth(CurrEvent=="StimShowing")')';
elDiff= diff(elevation(CurrEvent=="StimShowing")')';
nansum(abs(AzDiff))/nansum(abs(elDiff))

AzDiffNav= diff(azimuth(CurrEvent=="inTrial")')';
elDiffNav= diff(elevation(CurrEvent=="inTrial")')';
nansum(abs(AzDiffNav))/nansum(abs(elDiffNav))

% s=pwd;
% headposFactors= [mean(AzDiff),mean(elDiff),mean(AzDiffNav),mean(elDiffNav)];
% sTable = array2table(headposFactors);
% name=string(s)+'\SummaryStats\ReachheadposFactorsRad.csv';
% writetable(sTable,name);


%% Epoch eye movement plot
GazePos(GazePos==-1)=NaN;
GazeDiff= diff(GazePos')
GazeDiff=GazeDiff'
numOfTrial=50;
trialNum=Controller3DHeadContin(4,:);
trialNum=cellfun(@str2double,trialNum);

trialTime=Controller3DHeadContin(3,:);
trialTime=cellfun(@str2double,trialTime);
TrialEvents=CurrEvent(:,trialNum==numOfTrial);
MyIndex = find(contains(TrialEvents,'ObjDetected'));

moveTable=NaN(151,5000);
timeTable=NaN(151,5000);

figure()
for ii= 1:max(trialNum-1)
    numOfTrial=ii;
    Currgaze= GazeDiff(:,trialNum==numOfTrial);
    ThistrialTime= trialTime(:,trialNum==numOfTrial);
    TrialEvents=CurrEvent(:,trialNum==numOfTrial);
    MyIndex = find(contains(TrialEvents,'ObjDetected')); 
    ThistrialTime=ThistrialTime-ThistrialTime(MyIndex);
    
    M = movmean(sqrt(Currgaze(1,:).^2+Currgaze(2,:).^2+Currgaze(3,:).^2),25);
    ThistrialTime(end)=[];
    M(end)=[];
    plot(ThistrialTime,M,'k');
    hold on
    xline(ThistrialTime(MyIndex),'r');
    %align to index
    moveTable(ii,1:length(M))=M;
    timeTable(ii,1:length(ThistrialTime))=ThistrialTime;
end

%% Plot Average Eye Velocity
xlim([-5,5]);
moveTable(isnan(moveTable))=0;
times= linspace(-5,10,151);
step= (times(2)-times(1))/2;
ourMeans=[]
ourStes=[]
for ii= 1:length(times)
    thisArray =  moveTable(timeTable<times(ii)+step & timeTable>times(ii)-step);
    meanMove= nanmean(thisArray);
    ste= nanstd(thisArray)/sqrt(length(thisArray(~isnan(thisArray))));
    ourMeans=[ourMeans,meanMove];
    ourStes=[ourStes,ste];
end

figure()
plot(times,ourMeans.*60.*57.3,'color',myBlue, 'linewidth',2)
hold on
plot(times(times>0 & times<9.6),ourMeans(times>0 & times<9.6)*60*57.3,'color',myCyan, 'linewidth',2)
xline(0,'--','color',myRed, 'linewidth',2)
title('Eye movement in trial epochs','fontweight','bold','FontSize',TfontSize)
xlabel('Time(s)','fontweight','bold','FontSize',LfontSize)
ylabel('Average Eye Movement(degree/s)','fontweight','bold','FontSize',LfontSize)
width=500;
height=400;
x0=100;
y0=100;
set(gcf,'position',[x0,y0,width,height])
xticklabels({'Start','Detect','Locate'})
ax = gca;
xlim([-5,8]);
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
s=pwd;
name=string(s)+'\SummaryStats\ReacEyeMovementInTrialEpochs.jpg';
saveas(gcf,name)
name=string(s)+'\SummaryStats\ReacEyeMovementInTrialEpochs.svg';
saveas(gcf,name)


sTable = array2table(ourMeans);
name=string(s)+'\SummaryStats\ReacEyeVelocityEpochs.csv';
writetable(sTable,name);
%% Bar Graph EyeMovement during different types of navigation

thisTrialHitpoinst= cellfun(@str2double,trialHitpoints);
thisTrialeyeOpenness= cellfun(@str2double,eyeOpenness);
thisTrialevent=CurrEvent;

movment= Controller3DHeadContin(6:8,:);
movmentNew=cellfun(@str2double,movment);
DeltaMovement= movmentNew
mdiff = diff(movmentNew');
mdiff=mdiff';
EyeDiff= diff(thisTrialHitpoinst');
EyeDiff=EyeDiff';

%filter for navigate
thisTrialevent(:,1)=[];
mdiffNavigate= mdiff(:,thisTrialevent=="inTrial");
EyeDiffNavigate= EyeDiff(:,thisTrialevent=="inTrial");

this= mdiffNavigate(2,:)>.0001 & (mdiffNavigate(1,:)+mdiffNavigate(3,:))<.0001;
EyeVertNav= EyeDiffNavigate(:, this);
this= mdiffNavigate(2,:)<.0001 & (mdiffNavigate(1,:)+mdiffNavigate(3,:))>.0001;
EyeHorzNav= EyeDiffNavigate(:, this);

MeanVertNavEye= mean(abs(EyeVertNav'));
STDVertNavEye= std(abs(EyeVertNav'));

x=[1,2,3];
STEVertNavEye=STDVertNavEye/sqrt(length(EyeVertNav))

figure()
subplot(1,2,1)
bar(x,[MeanVertNavEye(1),MeanVertNavEye(3),MeanVertNavEye(2)]*100)     
xticklabels({'X','Y','Z'})
hold on
errorbar(x,[MeanVertNavEye(1),MeanVertNavEye(3),MeanVertNavEye(2)]*100, STEVertNavEye*100,'.','Color','k'); 
ylim([0,20])
grid on
title('Vertical Nav')
xlabel('Dimension')
ylabel('Eye Movement')

MeanHorzNavEye= mean(abs(EyeHorzNav'));
STDHorzNavEye= std(abs(EyeHorzNav'));
STEHorzNavEye=STDHorzNavEye/sqrt(length(EyeHorzNav));

subplot(1,2,2)
bar(x,[MeanHorzNavEye(1),MeanHorzNavEye(3),MeanHorzNavEye(2)]*100);     
xticklabels({'X','Y','Z'});
hold on
errorbar(x,[MeanHorzNavEye(1),MeanHorzNavEye(3),MeanHorzNavEye(2)]*100, STEHorzNavEye*100,'.','Color','k'); 
ylim([0,20])
grid on
title('Horizontal Nav')
xlabel('Dimension')
ylabel('Eye Movement')

s = pwd;
name=string(s)+'\SummaryStats\ReachNavEyeData.jpg';
saveas(gcf,name)

% Statistical analysis 
[h,p,ci,stats] =ttest2(EyeVertNav(1,:),EyeVertNav(2,:));
[h,p,ci,stats] =ttest2(EyeVertNav(1,:),EyeVertNav(3,:));
[h,p,ci,stats] =ttest2(EyeVertNav(2,:),EyeVertNav(3,:));

[h,p,ci,stats] =ttest2(EyeHorzNav(1,:),EyeHorzNav(2,:));
[h,p,ci,stats] =ttest2(EyeHorzNav(1,:),EyeHorzNav(3,:));
[h,p,ci,stats] =ttest2(EyeHorzNav(2,:),EyeHorzNav(3,:));

%Format and save data
NavEyeData(1,:)=MeanVertNavEye;
NavEyeData( 1, end+1:end+3 )=MeanHorzNavEye;

colNames = {'MeanVertNavEyeX','MeanVertNavEyeZ','MeanVertNavEyeY','MeanHorzNavEyeX','MeanHorzNavEyeZ','MeanHorzNavEyeY'};
rowNames = {'1'};
sTable = array2table(NavEyeData,'RowNames',rowNames,'VariableNames',colNames);
name=string(s)+'\SummaryStats\ReachNavEyeData.csv';
writetable(sTable,name);


%% Scatterplot of eye in head positions

GazePos=Controller3DHeadContin(12:14,:);
GazePos=cellfun(@str2double,GazePos);
GazePos(GazePos==-1)=NaN;
DotSize=20;
TfontSize=16;
LfontSize=12;
trialTime=cellfun(@str2double,Controller3DHeadContin(3,:));

figure()
x0=400;
y0=200;
width=1000;
height=400;
set(gcf,'position',[x0,y0,width,height])
subplot(1,2,1)
scatter(GazePos(1,CurrEvent=="inTrial"),GazePos(2,CurrEvent=="inTrial"),'.','SizeData',DotSize,'MarkerEdgeColor','b')
axis equal
xlabel('X','fontweight','bold','FontSize',LfontSize)
ylabel('Y','fontweight','bold','FontSize',LfontSize)
title('Eye Position During Navigation','fontweight','bold','FontSize',TfontSize)
grid on

ylim([-1,1]);
xlim([-1,1]);
xticks(-1:0.25:1);
yticks(-1:0.25:1);
subplot(1,2,2)
scatter(GazePos(1,CurrEvent=="StimShowing"),GazePos(2,CurrEvent=="StimShowing"),'.','SizeData',DotSize,'MarkerEdgeColor','r')
grid on
LfontSize= 12;
axis equal
xlabel('X','fontweight','bold','FontSize',LfontSize)
ylabel('Y','fontweight','bold','FontSize',LfontSize)
title('Eye Position During Search','fontweight','bold','FontSize',TfontSize)
ylim([-1,1]);
xlim([-1,1]);
xticks(-1:0.25:1);
yticks(-1:0.25:1);

s = pwd;
name=string(s)+'\SummaryStats\ReachGazePos.jpg';
saveas(gcf,name)

% save covariacnce matrix of the data
CovGazePosNav=nancov(GazePos(:,CurrEvent=="inTrial")');
colNames = {'X','Y','Z'};
sTable = array2table(CovGazePosNav,'VariableNames',colNames);
name=string(s)+'\SummaryStats\ReachCovGazePosNav.csv';
writetable(sTable,name);

CovGazePosSearch=nancov(GazePos(:,CurrEvent=="StimShowing")');
colNames = {'X','Y','Z'};
sTable = array2table(CovGazePosSearch,'VariableNames',colNames);
name=string(s)+'\SummaryStats\ReachCovGazePosSearch.csv';
writetable(sTable,name);

navEyePos=(GazePos(:,CurrEvent=="inTrial")');
colNames = {'X','Y','Z'};
sTable = array2table(navEyePos,'VariableNames',colNames);
name=string(s)+'\SummaryStats\ReachnavEyePos.csv';
writetable(sTable,name);
searchEyePos=(GazePos(:,CurrEvent=="StimShowing")');
colNames = {'X','Y','Z'};
sTable = array2table(searchEyePos,'VariableNames',colNames);
name=string(s)+'\SummaryStats\ReachsearchEyePos.csv';
writetable(sTable,name);
%% Plot Eye in head Heatmaps
pts = linspace(-1, 1, 41);
N = histcounts2(GazePos(2,CurrEvent=="inTrial"), GazePos(1,CurrEvent=="inTrial"), pts, pts);

figure()
x0=400;
y0=200;
width=1000;
height=400;
set(gcf,'position',[x0,y0,width,height])
subplot(1,2,1)
imagesc(pts, pts, N);
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');
% scatter(GazePos(1,CurrEvent=="inTrial"),GazePos(2,CurrEvent=="inTrial"),'.','SizeData',DotSize,'MarkerEdgeColor','b')
axis equal
xlabel('X','fontweight','bold','FontSize',LfontSize)
ylabel('Y','fontweight','bold','FontSize',LfontSize)
title('Eye Position During Navigation','fontweight','bold','FontSize',TfontSize)
grid on

ylim([-1,1])
xlim([-1,1])
xticks(-1:0.25:1);
yticks(-1:0.25:1);

subplot(1,2,2)
pts = linspace(-1, 1, 41);
N = histcounts2(GazePos(2,CurrEvent=="StimShowing"), GazePos(1,CurrEvent=="StimShowing"), pts, pts);
imagesc(pts, pts, N);
axis equal;
set(gca, 'XLim', pts([1 end]), 'YLim', pts([1 end]), 'YDir', 'normal');

LfontSize= 12;
axis equal
xlabel('X','fontweight','bold','FontSize',LfontSize)
ylabel('Y','fontweight','bold','FontSize',LfontSize)
title('Eye Position During Search','fontweight','bold','FontSize',TfontSize)
ylim([-1,1])
xlim([-1,1])
xticks(-1:0.25:1);
yticks(-1:0.25:1);
name=string(s)+'\SummaryStats\ReachGazePosHeatmap.jpg';
saveas(gcf,name)
%% Gaze Movement Scatterplot in epochs

GazePos(GazePos==-1)=NaN;
GazeDiff= diff(GazePos');
GazeDiff=GazeDiff';

figure()
x0=400;
y0=200;
width=1000;
height=400;

set(gcf,'position',[x0,y0,width,height])
subplot(1,2,1)
%CurrEvent(1)=[]
scatter(GazeDiff(1,CurrEvent=="inTrial"),GazeDiff(2,CurrEvent=="inTrial"),'.','SizeData',DotSize,'MarkerEdgeColor','b')
axis equal
xlabel('X','fontweight','bold','FontSize',LfontSize)
ylabel('Y','fontweight','bold','FontSize',LfontSize)
title('Eye movement During Navigation','fontweight','bold','FontSize',TfontSize)
grid on
ylim([-1,1])
xlim([-1,1])
xticks(-1:0.25:1);
yticks(-1:0.25:1);

subplot(1,2,2)
scatter(GazeDiff(1,CurrEvent=="StimShowing"),GazeDiff(2,CurrEvent=="StimShowing"),'.','SizeData',DotSize,'MarkerEdgeColor','r')
grid on
LfontSize= 12;
axis equal
xlabel('X','fontweight','bold','FontSize',LfontSize)
ylabel('Y','fontweight','bold','FontSize',LfontSize)
title('Eye movement During Search','fontweight','bold','FontSize',TfontSize)
ylim([-1,1])
xlim([-1,1])
xticks(-1:0.25:1);
yticks(-1:0.25:1);

s = pwd;
name=string(s)+'\SummaryStats\ReachGazeMovement.jpg';
saveas(gcf,name)

CovGazeMoveNav=nancov(GazeDiff(:,CurrEvent=="inTrial")');
colNames = {'X','Y','Z'};
sTable = array2table(CovnavGazeMoveNav,'VariableNames',colNames);
name=string(s)+'\SummaryStats\ReachCovGazeMoveNav.csv';
writetable(sTable,name);

CovGazeMoveSearch=nancov(GazeDiff(:,CurrEvent=="StimShowing")');
colNames = {'X','Y','Z'};
sTable = array2table(CovGazeMoveSearch,'VariableNames',colNames);
name=string(s)+'\SummaryStats\ReachCovGazeMoveSearch.csv';
writetable(sTable,name);
%% Saving it

navEyeMove=(GazeDiff(:,CurrEvent=="inTrial")');
colNames = {'X','Y','Z'};
nanvar(navEyeMove)
sTable = array2table(navEyeMove,'VariableNames',colNames);
name=string(s)+'\SummaryStats\navEyeMove.csv';
writetable(sTable,name);
searchEyeMove=(GazeDiff(:,CurrEvent=="StimShowing")');
colNames = {'X','Y','Z'};
sTable = array2table(searchEyeMove,'VariableNames',colNames);
name=string(s)+'\SummaryStats\searchEyeMove.csv';
writetable(sTable,name);

%%  Epoch head movement plot

Movement=Controller3DHeadContin(6:8,:);
Movement=cellfun(@str2double,Movement);
Movement= diff(Movement');
Movement=Movement';
trialNum=Controller3DHeadContin(4,:);
trialNum=cellfun(@str2double,trialNum);
trialTime=Controller3DHeadContin(3,:);
trialTime=cellfun(@str2double,trialTime);
moveTable=NaN(151,5000);
timeTable=NaN(151,5000);

figure()
for ii= 1:max(trialNum)-1
    numOfTrial=ii;
    Currmove= Movement(:,trialNum==numOfTrial);
    ThistrialTime= trialTime(:,trialNum==numOfTrial);
    TrialEvents=CurrEvent(:,trialNum==numOfTrial);
    MyIndex = find(contains(TrialEvents,'ObjDetected')); 
    ThistrialTime=ThistrialTime-ThistrialTime(MyIndex);
    M = movmean(sqrt(Currmove(1,:).^2+Currmove(2,:).^2+Currmove(3,:).^2),20);
    ThistrialTime(end)=[];
    M(end)=[];
    plot(ThistrialTime,M,'k');
    hold on
    xline(ThistrialTime(MyIndex),'r');
    %align to index
    moveTable(ii,1:length(M))=M;
    timeTable(ii,1:length(ThistrialTime))=ThistrialTime;
end

xlim([-5,10]);

moveTable(isnan(moveTable))=0;
times= linspace(-5,10,151);
step= (times(2)-times(1))/2;
ourMeans=[];
ourStes=[];
lengths=[];
for ii= 1:length(times)
    thisArray =  moveTable(timeTable<times(ii)+step & timeTable>times(ii)-step);
    for tt= 1:(30-sqrt(length(thisArray(~isnan(thisArray)))))
        thisArray=[thisArray;0]
    end
    meanMove= mean(thisArray);
    ste= nanstd(thisArray)/sqrt(length(thisArray(~isnan(thisArray))));
    ourMeans=[ourMeans,meanMove];
    ourStes=[ourStes,ste];
    lengths=[lengths,sqrt(length(thisArray(~isnan(thisArray))))];
end

figure()
plot(times,ourMeans*100,'color',myBlue, 'linewidth',2)
hold on
plot(times(times>0 & times<9.6),ourMeans(times>0 & times<9.6)*100,'color',myCyan, 'linewidth',2)
xline(0,'--','color',myRed, 'linewidth',2)
title('Movement in trial epochs','fontweight','bold','FontSize',TfontSize)
xlabel('Time(s)','fontweight','bold','FontSize',LfontSize)
ylabel('Average Velocity(m/s)','fontweight','bold','FontSize',LfontSize)
width=500;
height=400;
x0=100;
y0=100;
set(gcf,'position',[x0,y0,width,height])
xticklabels({'Start','Detect','Locate'})
ax = gca;
ylim([-0,.8]);
xlim([-5,9]);
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
name=string(s)+'\SummaryStats\ReachMovementInTrialEpochs.jpg';
saveas(gcf,name)
name=string(s)+'\SummaryStats\ReachMovementInTrialEpochs.svg';
saveas(gcf,name)


%%  Epoch Rotation plot

Movement=Controller3DHeadContin(9:11,:);
Movement=cellfun(@str2double,Movement);

Movement= diff(Movement');
Movement=Movement';

trialNum=Controller3DHeadContin(4,:);
trialNum=cellfun(@str2double,trialNum);

trialTime=Controller3DHeadContin(3,:);
trialTime=cellfun(@str2double,trialTime);
TrialEvents=CurrEvent(:,trialNum==numOfTrial);
MyIndex = find(contains(TrialEvents,'ObjDetected'));

moveTable=NaN(151,5000);
timeTable=NaN(151,5000);


figure()
for ii= 1:max(trialNum)-1
    numOfTrial=ii;
    Currmove= Movement(:,trialNum==numOfTrial);
    ThistrialTime= trialTime(:,trialNum==numOfTrial);
    TrialEvents=CurrEvent(:,trialNum==numOfTrial);
    MyIndex = find(contains(TrialEvents,'ObjDetected')); 
    ThistrialTime=ThistrialTime-ThistrialTime(MyIndex);
    M = movmean(sqrt(Currmove(1,:).^2+Currmove(2,:).^2+Currmove(3,:).^2),25);
    ThistrialTime(end)=[];
    M(end)=[];
    plot(ThistrialTime,M,'k');
    hold on
    xline(ThistrialTime(MyIndex),'r');
    %align to index
    moveTable(ii,1:length(M))=M;
    timeTable(ii,1:length(ThistrialTime))=ThistrialTime;
end

xlim([-5,10]);

moveTable(isnan(moveTable))=0;
times= linspace(-5,10,151);
step= (times(2)-times(1))/2;
ourMeans=[]
ourStes=[]
lengths=[]
for ii= 1:length(times)
    thisArray =  moveTable(timeTable<times(ii)+step & timeTable>times(ii)-step);
    for tt= 1:(30-sqrt(length(thisArray(~isnan(thisArray)))))
        thisArray=[thisArray;0]
    end
    meanMove= mean(thisArray);
    ste= nanstd(thisArray)/sqrt(length(thisArray(~isnan(thisArray))));
    ourMeans=[ourMeans,meanMove];
    ourStes=[ourStes,ste];
    lengths=[lengths,sqrt(length(thisArray(~isnan(thisArray))))];
end

figure()
plot(times,ourMeans*100*57,'color',myBlue, 'linewidth',2)
hold on
plot(times(times>0 & times<9.6),ourMeans(times>0 & times<9.6)*100*57,'color',myCyan, 'linewidth',2)
xline(0,'--','color',myRed, 'linewidth',2)
title('Rotation in trial epochs','fontweight','bold','FontSize',TfontSize)
xlabel('Time(s)','fontweight','bold','FontSize',LfontSize)
ylabel('AVG Head Rotation Speed(degree/s)','fontweight','bold','FontSize',LfontSize)
width=500;
height=400;
x0=100;
y0=100;
set(gcf,'position',[x0,y0,width,height])
xticklabels({'Start','Detect','Locate'})
ax = gca;
ylim([-0,90]);
xlim([-5,8]);
ax.XAxis.FontSize = ticksizeX;
ax.YAxis.FontSize = ticksizeY;
s=pwd;
name=string(s)+'\SummaryStats\ReachHeadRotInTrialEpochs.jpg';
saveas(gcf,name)
name=string(s)+'\SummaryStats\ReachHeadRotInTrialEpochs.svg';
saveas(gcf,name)

%%
Movement(:,1:6)=[]

figure()
subplot(3,1,1)
plot(Movement(1,:))
subplot(3,1,2)
plot(Movement(3,:))
subplot(3,1,3)
plot(Movement(2,:))

pdf('Poisson',Movement(1,:));
[f,xi] = ksdensity(Movement(1,:));
figure()
subplot(3,1,1)
plot(xi,f);
subplot(3,1,2)
[f,xi] = ksdensity(Movement(3,:));
plot(xi,f);
subplot(3,1,3)
[f,xi] = ksdensity(Movement(2,:));
plot(xi,f);




%% Head rotation position Trial

rot=Controller3DHeadContin(9:11,:);
rot=cellfun(@str2double,rot);
trialNum=Controller3DHeadContin(4,:);
trialNum=cellfun(@str2double,trialNum);

figure()
x0=400;
y0=200;
width=1000;
height=400;
DotSize=5;
set(gcf,'position',[x0,y0,width,height])
subplot(1,2,1)
scatter3(rot(1,CurrEvent=="inTrial"),rot(3,CurrEvent=="inTrial"),rot(2,CurrEvent=="inTrial"),'.','SizeData',DotSize,'MarkerEdgeColor',myBlue)
axis equal
xlabel('X','fontweight','bold','FontSize',LfontSize);
ylabel('Y','fontweight','bold','FontSize',LfontSize);
title('Head rotation During Navigation','fontweight','bold','FontSize',TfontSize)
grid on
ylim([-1,1])
xlim([-1,1])
zlim([-1,1])
xticks(-1:0.5:1);
yticks(-1:0.5:1);
zticks(-1:0.5:1);

subplot(1,2,2)
scatter3(rot(1,CurrEvent=="StimShowing"),rot(3,CurrEvent=="StimShowing"),rot(2,CurrEvent=="StimShowing"),'.','SizeData',DotSize,'MarkerEdgeColor',myRed)
%scatter3(GazePos(1,CurrEvent=="StimShowing"),GazePos(2,CurrEvent=="StimShowing"),GazePos(3,CurrEvent=="StimShowing"),20,trialTime(CurrEvent=="StimShowing"), 'filled','o')
grid on
LfontSize= 12;
axis equal
xlabel('X','fontweight','bold','FontSize',LfontSize);
ylabel('Y','fontweight','bold','FontSize',LfontSize);
title('Head rotation During Search','fontweight','bold','FontSize',TfontSize)
ylim([-1,1])
xlim([-1,1])
zlim([-1,1])
xticks(-1:0.5:1);
yticks(-1:0.5:1);
zticks(-1:0.5:1);
% Save it
s = pwd;
name=string(s)+'\SummaryStats\ReachHeadrotPos.jpg';
saveas(gcf,name);
s = pwd;
name=string(s)+'\SummaryStats\ReachHeadrotPos.svg';
saveas(gcf,name);

%%  Extract and save vertical and horizontal factors of navigation in epochs
rot=Controller3DHeadContin(9:11,:);
rot=cellfun(@str2double,rot);
trialNum=Controller3DHeadContin(4,:);
trialNum=cellfun(@str2double,trialNum);
horzFactorNav= real(mean(sqrt((rot(1,CurrEvent=="StimShowing")+rot(3,CurrEvent=="StimShowing")))));
vertFactorNav= real(mean(sqrt((rot(1,CurrEvent=="StimShowing")+rot(2,CurrEvent=="StimShowing")))));

horzFactorSearch= real(mean(sqrt((rot(1,CurrEvent=="inTrial")+rot(3,CurrEvent=="inTrial")))));
vertFactorSearch= real(mean(sqrt((rot(1,CurrEvent=="inTrial")+rot(2,CurrEvent=="inTrial")))));

headposFactors= [horzFactorNav,vertFactorNav,horzFactorSearch,vertFactorSearch];
colNames = {'horzFactorNav','vertFactorNav','horzFactorSearch','vertFactorSearch'};
sTable = array2table(headposFactors,'VariableNames',colNames);
name=string(s)+'\SummaryStats\headposFactors.csv';
writetable(sTable,name);

s=pwd;
headposSearch= [rot(1,CurrEvent=="StimShowing");rot(3,CurrEvent=="StimShowing");rot(2,CurrEvent=="StimShowing")];
sTable = array2table(headposSearch);
name=string(s)+'\SummaryStats\headposSearch.csv';
writetable(sTable,name);

headposNav= [rot(1,CurrEvent=="inTrial");rot(3,CurrEvent=="inTrial");rot(2,CurrEvent=="inTrial")];
sTable = array2table(headposNav);
name=string(s)+'\SummaryStats\headposNav.csv';
writetable(sTable,name);


%% Factors during reach

[azimuth,elevation,r] =  cart2sph(rot(1,:),rot(3,:),rot(2,:));

nanmean(azimuth)*90;
nanmean(elevation)*90;

AzDiff= diff(azimuth(CurrEvent=="StimShowing")')';
elDiff= diff(elevation(CurrEvent=="StimShowing")')';
nansum(abs(AzDiff))/nansum(abs(elDiff))

AzDiffNav= diff(azimuth(CurrEvent=="inTrial")')';
elDiffNav= diff(elevation(CurrEvent=="inTrial")')';
nansum(abs(AzDiffNav))/nansum(abs(elDiffNav));

headposFactors= [AzDiff,elDiff,AzDiffNav,elDiffNav];
sTable = array2table(headposFactors);
name=string(s)+'\SummaryStats\ReachheadposFactorsRad.csv';
writetable(sTable,name);

%% Head rotation movement during Epochs in polar coordinates

[azimuth,elevation,r]= cart2sph(rot(1,:),rot(2,:),rot(3,:));
AzChange= (diff(azimuth'))';
ElChange= (diff(elevation'))';
r= (diff(r'))';
[thisX,thisY,thisZ]= sph2cart(AzChange,ElChange,r);

figure()
x0=400;
y0=200;
width=1000;
height=400;
DotSize=5;
set(gcf,'position',[x0,y0,width,height])
subplot(1,2,1)
scatter(AzChange(CurrEvent=="inTrial"),ElChange(CurrEvent=="inTrial")'.','SizeData',DotSize,'MarkerEdgeColor','b')

axis equal
xlabel('X','fontweight','bold','FontSize',LfontSize)
ylabel('Y','fontweight','bold','FontSize',LfontSize)
title('During Navigation','fontweight','bold','FontSize',TfontSize)
grid on
ylim([-0.03,0.03])
xlim([-0.03,0.03])
zlim([-0.03,0.03])
xticks(-1:0.02:1);
yticks(-1:0.02:1);
zticks(-1:0.02:1);
subplot(1,2,2)
scatter(AzChange(CurrEvent=="StimShowing"),ElChange(CurrEvent=="StimShowing")'.','SizeData',DotSize,'MarkerEdgeColor','b')

grid on
LfontSize= 12;
axis equal
xlabel('X','fontweight','bold','FontSize',LfontSize)
ylabel('Y','fontweight','bold','FontSize',LfontSize)
title('During Search','fontweight','bold','FontSize',TfontSize)
ylim([-0.03,0.03])
xlim([-0.03,0.03])
zlim([-0.03,0.03])
xticks(-1:0.02:1);
yticks(-1:0.02:1);
zticks(-1:0.02:1);

s = pwd;
name=string(s)+'\SummaryStats\ReachHeadrotMove.jpg';
saveas(gcf,name)

CovheadMoveNav=nancov(rotChange(:,CurrEvent=="inTrial")');
colNames = {'X','Y','Z'};
sTable = array2table(CovheadMoveNav,'VariableNames',colNames);
name=string(s)+'\SummaryStats\ReachCovheadMoveNav.csv';
writetable(sTable,name);

CovheadMoveSearch=nancov(rotChange(:,CurrEvent=="StimShowing")');
colNames = {'X','Y','Z'};
sTable = array2table(CovheadMoveSearch,'VariableNames',colNames);
name=string(s)+'\SummaryStats\ReachCovheadMoveSearch.csv';
writetable(sTable,name);

%% Head rotation position + eye rotation position

%make matrix same size as rot
b= ones(size(rot));
r= []
r = vrrotvec(rot(:,2),b(:,2))

for ii=2:length(rot)
     r(ii,:)=vrrotvec(rot(:,ii),b(:,ii));
end

ThisGazePos = GazePos(flip([1 3 2]), :);

figure()
x0=400;
y0=200;
width=1000;
height=400;
DotSize=5;
set(gcf,'position',[x0,y0,width,height])
subplot(1,2,1)
scatter3(rot(1,CurrEvent=="inTrial"),rot(3,CurrEvent=="inTrial"),rot(2,CurrEvent=="inTrial"),'.','SizeData',DotSize,'MarkerEdgeColor','black')
for ii=1:10000
    hold on
    if CurrEvent(ii)=="inTrial"
        m = vrrotvec2mat(r(ii,:));
        newGaze=m*ThisGazePos(:,ii);
        scatter3(rot(1,ii)+ newGaze(1),rot(3,ii)+newGaze(2),rot(2,ii)+newGaze(3),'.','SizeData',DotSize,'MarkerEdgeColor','b');
    end
end

axis equal
xlabel('X','fontweight','bold','FontSize',LfontSize)
ylabel('Y','fontweight','bold','FontSize',LfontSize)
title('head rotation During Navigation','fontweight','bold','FontSize',TfontSize)
grid on
ylim([-1,1])
xlim([-1,1])

subplot(1,2,2)
scatter3(rot(1,CurrEvent=="StimShowing"),rot(3,CurrEvent=="StimShowing"),rot(2,CurrEvent=="StimShowing"),'.','SizeData',DotSize,'MarkerEdgeColor','r')
grid on
LfontSize= 12;
axis equal
xlabel('X','fontweight','bold','FontSize',LfontSize)
ylabel('Y','fontweight','bold','FontSize',LfontSize)
title('Head rotation During Search','fontweight','bold','FontSize',TfontSize)
ylim([-1,1])
xlim([-1,1])


%% PLot Gaze position in arena for trial

Controller3DHeadContin=(Controller3DHeadContin);
trialHitpoints=Controller3DHeadContin(20:22,:);
Controller3DHead= xlsread('Data-4-Controller3DHead.csv');

%pick trialNum
trialNumber=42;
TrialNum=Controller3DHeadContin(4,:);
CurrEvent=Controller3DHeadContin(5,:);
eyeOpenness=Controller3DHeadContin(38:39,:);
trialTime=cellfun(@str2double,Controller3DHeadContin(3,:));
trialNumberstring=string(trialNumber);
thisTrialHitpoinst=Controller3DHeadContin(20:22,TrialNum==trialNumberstring);
thisTrialevent=CurrEvent(:,TrialNum==trialNumberstring);
thisTrialeyeOpenness=eyeOpenness(:,TrialNum==trialNumberstring);

thisTrialHitpoinst= cellfun(@str2double,thisTrialHitpoinst);
thisTrialeyeOpenness= cellfun(@str2double,thisTrialeyeOpenness);
thisTrialeyeOpennessmean=mean(thisTrialeyeOpenness);

thisTrialHitpoinst=thisTrialHitpoinst(:,thisTrialeyeOpennessmean>0.98);
thisTrialevent=thisTrialevent(:,thisTrialeyeOpennessmean>0.98);
Controller3DHead(:,100:end)=[];
locAtDetect=Controller3DHead(16:18,:);
loctarget=Controller3DHead(2:4,:);
locDot=Controller3DHead(5:7,:);

%then plot it
figure()
plot3(thisTrialHitpoinst(1,thisTrialevent=="inTrial"),thisTrialHitpoinst(3,thisTrialevent=="inTrial"),thisTrialHitpoinst(2,thisTrialevent=="inTrial"));
hold on;
scatter3(thisTrialHitpoinst(1,thisTrialevent=="inTrial"),thisTrialHitpoinst(3,thisTrialevent=="inTrial"),thisTrialHitpoinst(2,thisTrialevent=="inTrial"),20,trialTime(thisTrialevent=="inTrial"), 'filled','o');
colormap(cool);
grid on;
plot3(1,1,1);
axis equal;
%plot target location
scatter3(locAtDetect(1,trialNumber),locAtDetect(3,trialNumber),locAtDetect(2,trialNumber),50,'filled','o','red');
scatter3(loctarget(1,trialNumber),loctarget(3,trialNumber),loctarget(2,trialNumber),50,'filled','o','blue');
scatter3(locDot(1,trialNumber),locDot(3,trialNumber),locDot(2,trialNumber),50,'filled','o','black');
%plot sub start point

xlim([-5 5])
ylim([-5 5])
zlim([-1 9])





