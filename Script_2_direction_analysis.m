%%
% clearvars
clearvars -except cdelx_all cdely_all sdelx_all sdely_all 
close all
clc

%You need to run script 1: track analysis before this script

%IMPORTANT:
%mannually input file name, xlsread, and column values

%INPUT file name of each experiment repeat (one control and one tnnt2a/sih)
%Make sure it's the same name used in previous script
ctrl = '170927_ctrl_analysis.xlsx';
sih = '170927_sih_analysis.xlsx';
delx_range = 'A:A'; %tip cells,'A:A'; vSIV cell,'C:C'; branch cell,'E:E'
dely_range = 'B:B'; %tip cells,'B:B'; vSIV cell,'D:D'; branch cell,'F:F'

%OUTPUT filenmae, export polarhistogram values of each repeat
Output_filename = 'tip_angular_counts_Norm_45Degree_every_repeat.xlsx';
%number of repeat, the first repeat =1; other repeats = any other number
repeat = 2; 
%First experiment repeat 'A:A, second 'B:B', third 'C:C'...etc
data_location = 'A:A';

sheet = 'Step deltax_deltay';

cdelx =xlsread(ctrl,sheet,delx_range); % read excel, ctrl delta x
cdely =xlsread(ctrl,sheet,dely_range); % read excel, ctrl delta y
sdelx =xlsread(sih,sheet,delx_range);  % read excel, sih delta x
sdely =xlsread(sih,sheet,dely_range);  % read excel, sih delta y
cdelx= rmmissing(cdelx); %remove missing value
cdely= rmmissing(cdely); 
sdelx= rmmissing(sdelx); 
sdely= rmmissing(sdely); 

%store delx & dely in different repeat together
if repeat == 1
    cdelx_all = [];
    cdely_all = [];
    sdelx_all = [];
    sdely_all = [];
    
    cdelx_all = [cdelx_all;cdelx];
    cdely_all = [cdely_all;cdely];
    sdelx_all = [sdelx_all;sdelx];
    sdely_all = [sdely_all;sdely];
else
    cdelx_all = [cdelx_all;cdelx];
    cdely_all = [cdely_all;cdely];
    sdelx_all = [sdelx_all;sdelx];
    sdely_all = [sdely_all;sdely];
end

%calculate theta
c = atan2(cdely,cdelx); %theta of control
s = atan2(sdely,sdelx);

% group angles
%the number of each angle is normalised to total number of angles
hc =polarhistogram(c,'BinWidth',pi/4,'Normalization','probability','DisplayStyle','stairs');
hold on;
%get relative number of each angle
c_bars = hc.Values

hs =polarhistogram(s,'BinWidth',pi/4,'Normalization','probability','DisplayStyle','stairs');
hold on;
%get relative number of each angle
s_bars = hs.Values 

ctrl_angle_Norm = c_bars'
sih_angle_Norm = s_bars'

writematrix(ctrl_angle_Norm,Output_filename,'Sheet',1,'Range',data_location)
writematrix(sih_angle_Norm,Output_filename,'Sheet',2,'Range',data_location)

%%
Output_filename = 'Tip_angular_counts_Norm_45Degree_every_repeat.xlsx';

header1= {'ctrl_delx','ctrl_dely','sih_delx','sih_dely'};
write_xls({cdelx_all,cdely_all,sdelx_all,sdely_all} ,Output_filename, 'Tip cell_all delta', header1);

%%
%calculate theta
control = atan2(cdely_all,cdelx_all); %theta of control (all repeats)
tnnt2a = atan2(sdely_all,sdelx_all);  %theta of sih (all repeats)

%OUTPUT figure name of scatter plot
fig_name_1 = 'ctrl_sih_pool_scatter_plot.svg'
fig_name_2 = 'ctrl_sih_pool_polarhitogram.svg'

%SCATTER PLOT
SP_title = 'SIV cells'
SP_xlabel = 'step delta x'
SP_ylabel = 'step delta y'

xCenter = 0;
yCenter = 0;
radius = 5;

fig1 = figure(100)
subplot(1,2,1);
scatter(cdelx_all,cdely_all,'k');
xlim([-25 25]);  %tip & vSIV [-15 15]; branch [-25 25]
ylim([-25 25]);
hold on;
% orange '[0.8500, 0.3250, 0.0980]'
viscircles([xCenter yCenter],radius,'color','m');
title(SP_title)
xlabel(SP_xlabel) 
ylabel(SP_ylabel) 
legend('control','Location','northwest','FontSize',40)
ax = gca;
ax.FontSize = 30;
grid on;
ax.GridAlpha = 0.5

subplot(1,2,2);
scatter(sdelx_all,sdely_all,'k');
xlim([-25 25]);
ylim([-25 25]);
hold on;
viscircles([xCenter yCenter],radius,'color','m');
title(SP_title)
xlabel(SP_xlabel) 
ylabel(SP_ylabel) 
legend('tnnt2a','Location','northwest','FontSize',40)
ax = gca;
ax.FontSize = 30;
grid on;
ax.GridAlpha = 0.5
% saveas(fig1,fig_name_1)

%PLOT POLAR HISTOGRAM
fig2 = figure(200)
hc=polarhistogram(control,'BinWidth',pi/8,'Normalization','probability','DisplayStyle','stairs');
hold on;

hs=polarhistogram(tnnt2a,'BinWidth',pi/8,'Normalization','probability','DisplayStyle','stairs');
hold on;

hc.LineWidth = 5;
hc.EdgeColor = 'k';
hs.LineWidth = 5;
hs.EdgeColor = 'm';

% title('Tip cells')
legend({'control','tnnt2a'},'Location','northwest')
ax = gca;
ax.FontSize = 30;
grid on;
ax.GridAlpha = 1
thetaticks(0:45:315)

% saveas(fig2,fig_name_2)

counts1 = hc.Values
counts2 = hs.Values


