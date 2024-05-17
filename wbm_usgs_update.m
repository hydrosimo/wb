%% WATER BALANCE MODEL THORNTWAITE (USGS)

clear all
close all
warning('off','all')

%% Loading data and variables (1987-2020)

% Upstream basin
load ETr_enza.dat  % mm (Hargreaves formula)
prec = readmatrix('p_monte.csv');
tmean = flipud(readmatrix('temp_diurne.csv'));  
runoff_obs = flipud(readmatrix('runoff_obs.csv'));  % observed data (Annali)
runoff_obs(isnan(runoff_obs))=0;
train = 3.3;       % °C, McCabe and Wolock, 1999
tsnow = -1;        % °C, > 1000 m, USGS
drofrac = 0.05;    % -, fraction of Prain that becomes runoff directly, USGS
meltmax = 0.5;     % -, maximum melting, McCabe and Wolock, 1999; Wolock and McCabe,1999
STC = 150;         % mm, soil-moisture capacity, USGS
rfactor = 0.5;     % -, run-off factor, Wolock and McCabe, 1999

%% Water-Balance programming

% snow
psnow = prec.*((train-tmean)/(train-tsnow));
psnow(psnow<0) = 0;  % snostor, mm
prain = prec-psnow;
prain(prain<0) = 0;

% runoff
dro = prain*drofrac; % direct runoff, mm
prem = prain-dro;    % precipitation remained, mm

% snow melting
snowstorage = (((tmean-tsnow)/(train-tsnow))*meltmax);  % -, snow melt fraction
snowstorage(snowstorage<0) = 0;
snowstorage(snowstorage>meltmax) = meltmax;
snowmelt = psnow.*snowstorage;   % mm, amount of snow melted in a month
ptot = snowmelt+prem;            % mm, total liquid water input to the soil

% evapotranspiration, soil-moisture storage (Alley W.M.)
PET = ETr_enza;            % flipud() to invert matrix rows (Hargreaves)
for i = 1:length(ptot)     % ptot-PET difference
    for j = 1:size(ptot,2)
        ptotPET(i,j) = ptot(i,j)-PET(i,j);
    end
end
% -------------------------------------------------------------------------
if  ptot(1,1)>PET(1,1)     % mm, initial soil-moisture storage
    AET(1,1) = PET(1,1);
    S(1,1) = (ptot(1,1)-PET(1,1));
else
    WT(1,1) = (abs(ptot(1,1)-PET(1,1)))*(S(1,1)/STC);
    S(1,1) = ptot(1,1)*exp(-(PET(1,1)-ptot(1,1))/STC)-WT(1,1);
    AET(1,1) = min(ptot(1,1)+WT(1,1),PET(1,1));
end 
for j = 2:size(ptot,2)     % mm, first year
    if ptot(1,j)>PET(1,j) 
        AET(1,j) = PET(1,j);
        replenish(1,j) = ptotPET(1,j);   % replenish water (ptot-PET)
        surplus(1,j) = max(S(1,j-1)+replenish(1,j)-STC,0);
        S(1,j) = min(ptot(1,j)-PET(1,j),STC);
    else 
        WT(1,j) = (abs(ptot(1,j)-PET(1,j)))*(S(1,j-1)/STC);
        S(1,j) = S(1,j-1)-WT(1,j);
        AET(1,j) = min(ptot(1,j)+WT(1,j),PET(1,j));
    end
end
if ptot(2,1)>PET(2,1)      % mm, january second year
    AET(2,1) = PET(2,1);
    replenish(2,1) = ptotPET(2,1);
    surplus(2,1) = max(S(1,12)+replenish(2,1)-STC,0);
    S(2,1) = min(ptot(2,1)-PET(2,1),STC);
else
    WT(2,1) = (abs(ptot(2,1)-PET(2,1)))*(S(2,1)/STC);
    S(2,1) = S(1,12)-WT(2,1);
    AET(2,1) = min(ptot(2,1)+WT(2,1),PET(2,1));
end
for i = 2:length(ptot)-1   % mm, total simulation
    for j = 2:size(ptot,2)
        if ptot(i,j)>PET(i,j)
            AET(i,j) = PET(i,j);
            replenish(i,j) = ptotPET(i,j);
            surplus(i,j) = max(S(i,j-1)+replenish(i,j)-STC,0);
            S(i,j) = min(ptot(i,j)-PET(i,j),STC);
        else
            WT(i,j) = (abs(ptot(i,j)-PET(i,j)))*(S(i,j-1)/STC);
            S(i,j) = S(i,j-1)-WT(i,j);
            AET(i,j) = min(ptot(i,j)+WT(i,j),PET(i,j));
        end
        if ptot(i+1,1)>PET(i+1,1)
            replenish(i+1,1) = ptotPET(i+1,1);
            surplus(i+1,1) = max(S(i,12)+replenish(i+1,1)-STC,0);
            AET(i+1,1) = PET(i+1,1);
            S(i+1,1) = min(ptot(i+1,1)-PET(i+1,1),STC);
        else
            WT(i+1,1) = (abs(ptot(i+1,1)-PET(i+1,1)))*(S(i,12)/STC);
            S(i+1,1) = S(i,12)-WT(i+1,1);
            AET(i+1,1) = min(ptot(i+1,1)+WT(i+1,1),PET(i+1,1));
        end
    end
end
for j = 2:size(ptot,2)     % mm, last year
    if ptot(length(ptot),j)>PET(length(ptot),j)
        AET(length(ptot),j) = PET(length(ptot),j);
        replenish(length(ptot),j) = ptotPET(length(ptot),j);
        surplus(length(ptot),j) = max(S(length(ptot),j-1)+replenish(length(ptot),j)-STC,0);
        S(length(ptot),j) = min(ptot(length(ptot),j)-PET(length(ptot),j),STC);
    else
        WT(length(ptot),j) = (abs(ptot(length(ptot),j)-PET(length(ptot),j)))*(S(length(ptot),j-1)/STC);
        S(length(ptot),j) = S(length(ptot),j-1)-WT(length(ptot),j);
        AET(length(ptot),j) = min(ptot(length(ptot),j)+WT(length(ptot),j),PET(length(ptot),j));
    end
end
deficit = PET-AET;

% run-off 
RO(1,1) = rfactor*surplus(1,1);   % 50% of surplus becomes runoff
for j = 2: size(ptot,2)           % mm, first year
    if surplus(1,j-1) == 0
        RO(1,j) = rfactor*surplus(1,j);
    else
        RO(1,j) = rfactor*surplus(1,j)+RO(1,j-1);   % the remaining surplus is carried over to the following month 
    end
end
for i = 2:length(ptot)            % mm, total simulation
    for j = 2:size(ptot,2)
        if surplus(i-1,12) == 0
            RO(i,1) = rfactor*surplus(i,1);
        else
            RO(i,1) = rfactor*surplus(i,1)+RO(i-1,12);
        end
        if surplus(i,j-1) == 0
            RO(i,j) = rfactor*surplus(i,j);
        else
            RO(i,j) = rfactor*surplus(i,j)+RO(i,j-1);
        end
    end
end
ROtot = dro+RO;                   % mm, total runoff        
ROfinal = [ROtot sum(ROtot,2)];   % mm, total runoff
HR = (mean(ROtot))/1000;          % m, monthly runoff depth
HP = (mean(prec))/1000;           % m, monthly precipitation depth
HR_12_19 = (mean(ROtot(21:34,:)))/1000;          % m, monthly runoff depth
HP_12_19 = (mean(prec(21:34,:)))/1000;           % m, monthly precipitation depth

% runoff coefficient
phi = HR./HP;
phi_12_19 = HR_12_19./HP_12_19;
phi_obs = 0.61;                        % annali idrologici, pt. II,  p. 91, 2019 (2012-2019)
phi_sim = 1.0*min(phi,1);              % simulated (1987-2020)
phi_sim_12_19 = 1.0*min(phi_12_19,1);  % simulated (2012-2019)

%% Saving data:

runoff = fopen('runoff_sim.dat','wt');
fprintf(runoff,'%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n',ROfinal.');
fclose(runoff);

%% Plots: 

% Date settings
lw = 1.0;
minx = 1;  
maxx = size(ptot,2);
miny = -500;
maxy = 500;

% Figure settings
figure(1)
set(1,'Units','centimeters','Position',[1,1,35,19]);
set(1,'defaultaxesfontname', 'CMU Serif');
set(1,'defaulttextfontname', 'CMU Serif');
set(gca,'FontSize',20,'Box','on');
diffsmooth = smooth(mean(ROtot(21:34,:)-runoff_obs));
plot(diffsmooth,'r','LineWidth',4*lw), hold on
vxt = minx:1:maxx;
vyt = miny:100:maxy;
set(gca,'xtick',vxt,'ytick',vyt,'XMinorTick','off','YMinorTick','on');
axis([minx maxx miny maxy]);
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
xtickangle(90);
ylabel('mean runoff gap (mm)');
xlabel('months');
legend('observed and simulated mean runoff gap','Location','NorthEast','NumColumns',1);
legend('boxoff')
grid on
set(gca,'PlotBoxAspectRatio',[2.5 1.5 2],'FontSize',20,'Ticklength',[0.015 0.00],'LineWidth',1.5,'Layer','top');

orient(figure(1),'landscape')
saveas(1,'figone','svg')
system('inkscape --export-filename=FIG_wbm_ink.pdf figone.svg')
system('rm figone.svg')
system('pdfcrop FIG_wbm_ink.pdf FIG_wbm_cropped.pdf')
system('rm FIG_wbm_ink.pdf')

fprintf(1,'\n')
fprintf(1,'-----------------------------------------------------------------------\n')
fprintf(1,'MEAN ANNUAL RUNOFF COEFFICIENT (Enza drainage basin closed at Vetto):\n')
fprintf(1,'- observed monthly runoff depth,  from 2012 to 2019 = %5.2f\n',phi_obs)
fprintf(1,'- simulated monthly runoff depth, from 2012 to 2019 = %5.2f\n',mean(phi_sim_12_19))
fprintf(1,'- simulated monthly runoff depth, from 1987 to 2020 = %5.2f\n',mean(phi_sim))
fprintf(1,'-----------------------------------------------------------------------\n')
fprintf(1,'\n')

