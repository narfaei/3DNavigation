%% Three session Analysis
clear;
clc;
%% Define variables
%Define starting Probabilities
pHorizontal=0.3;
pUp=0.3;
pDown=0.3;
%Define Gamma
Gamma=0.8;

%% get the data
Controller3DHead= xlsread('Data-4-Controller3DHead.csv');
%% clean format Data
Controller3DHead(:,101:end)=[];
Controller3DHead=Controller3DHead';
Controller3DHeadError=abs(Controller3DHead(:,2:4) - Controller3DHead(:,5:7));
CartController3DHeadError=(Controller3DHead(:,2:4) - Controller3DHead(:,5:7));
[CleanedController3DHeadError,Controller3Dindecies] = rmoutliers(Controller3DHeadError,'mean');
Controller3DHeadError=mean(CleanedController3DHeadError);
Controller3DHead(Controller3Dindecies,:)=[];
Controller3DHeadSTE=Controller3DHead(:,2:4) - Controller3DHead(:,5:7);
Controller3DHeadSTE= std(Controller3DHeadSTE) / sqrt(length(Controller3DHeadSTE))
Controller3DHeadSTE = Controller3DHeadSTE(:,[1 3 2]);

x=[1,2,3];

%% 1-----  Controller Errors compared to the actual location of stimuli (stim location as [0,0,0])
xdiffCentered(:,1)=Controller3DHead(:,5)-Controller3DHead(:,2);
ydiffCentered(:,1)=Controller3DHead(:,7)-Controller3DHead(:,4);
zdiffCentered(:,1)=Controller3DHead(:,6)-Controller3DHead(:,3);

edges=-3:.2:3; % Define Edges 20cm

% make a 25 by 25 3 Dim Grid and bin data in there
[~,~,binnedX]=histcounts(xdiffCentered,edges);
[~,~,binnedY]=histcounts(ydiffCentered,edges);
[~,~,binnedZ]=histcounts(zdiffCentered,edges);
combined=binnedX;
combined(:,2)=binnedY;
combined(:,3)=binnedZ;

DataHist(:,1)=histcounts(binnedX,1:1:30);
DataHist(:,2)=histcounts(binnedY,1:1:30);
DataHist(:,3)=histcounts(binnedZ,1:1:30);

Errors=DataHist;

%% search for minimum

options = optimset('Display','iter','PlotFcns',@optimplotfval);
options.MaxIter=100;
X0=[Gamma,pHorizontal,pUp,pDown];

handle=@SRFit3D;
[result,optError]=fminsearch(handle,X0,options,DataHist);

% Save
s = pwd;
colNames = {'Gamma','pHorizontal','pUp','pDown'};
thisresult = array2table(result,'VariableNames',colNames);

name=string(s)+'\SummeryStats\SRFitResultsNew.csv';
writetable(thisresult,name);

name=string(s)+'\SummeryStats\SRFitErrorNew.csv';
colNames = {'FitError'};

thisError = array2table(squeeze(optError),'VariableNames',colNames);
writetable(thisError,name);

%% Plot results to compare

%lets plot to check
figure();
subplot(3,2,1);
bar(Errors(:,1));
ylim([0 .2]);
grid on;
subplot(3,2,2);
bar(Model(:,1));
ylim([0 .2]);
grid on;

subplot(3,2,3);
bar(Errors(:,2));
ylim([0 .2]);
grid on;
subplot(3,2,4);
bar(Model(:,2));
ylim([0 .2]);
grid on;

subplot(3,2,5);
bar(Errors(:,3));
ylim([0 .2]);
grid on;
subplot(3,2,6);
bar(Model(:,3));
ylim([0 .2]);
grid on;



     