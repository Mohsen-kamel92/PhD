%% ===== Build the target torque traces from a maximal concentric contraction =====
%
% PURPOSE
%   Turn a participant's maximal concentric dorsiflexion into the torque traces
%   they follow during the experiment, at each target intensity and for both
%   contraction modes.
%
% INPUT
%   A Spike2 recording of one maximal concentric contraction, exported to .mat
%   with SaveMAT.s2s. Must contain Torque and Angle channels.
%
% OUTPUT
%   <chosen folder>/ConMatch_P1_<pct>.txt   or   EccMatch_P1_<pct>.txt
%   Two columns, time against torque, read by the Spike2 protocol into its XY
%   view and by Stage 01 of the analysis pipeline.
%
% METHOD
%   The rotation phase is located from the points where the joint angle stops
%   changing. The passive torque-angle relation is fitted and subtracted, the
%   active component is scaled to each target intensity, and the passive
%   component is added back, so the trace represents a torque the participant
%   can actually produce rather than a fraction of a total. The concentric trace
%   is reversed in time to give the eccentric trace, so both modes follow the
%   same torque at the same joint angle. Resting torque at the short and the
%   long position is appended at either end.
%
%   Only the rotation phase is scaled. Resting torque is left untouched, since
%   scaling it would misrepresent the passive load.
%
% ON THE FILE NAMES
%   P1 or P2 reflects how many scaling factors were entered in one run. It does
%   not identify the participant, which is why the analysis matches ramp files
%   on contraction mode and intensity alone.
%
% DEPENDENCIES
%   Signal Processing Toolbox   findchangepts, butter, filtfilt
%   Wcu, Auto_RestTorque, Normalize_Torque, dlg_SetTrace

clc
clear
close all

%% ---------------------------
%  Paths and input file
%% ---------------------------
% Helper functions live one level up, in the shared functions folder.
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'functions'));

%open MAT file of maximal concentric contraction from spike
[SpikeFile, pth] = uigetfile('*.mat');
load([pth SpikeFile]);

Frequency = 1/Torque.interval; %I'd use 200Hz would be enough

%% Find the abrut changes, index (1) where torque changes as depends on trigger threshold
% index(2) where Angle stop moving --> flexion movement
% index(3) where Angle stop moving again --> ext movement
Tor = Torque.values;
Ang = Angle.values;
%% Filtering torque so then plotting in spike is easy
CutOff = 7;
Filt_order = 4;

%correct order and cut off
[cu,order] = Wcu(Filt_order,CutOff,'low',Frequency);
[b,a] = butter(order, cu/(0.5*Frequency),'low');

Tor = filtfilt(b,a,Tor);
Ang = filtfilt(b,a,Ang);

%% Find changing pts based on angle

pts = findchangepts(diff(Ang),"MaxNumChanges",4,'Statistic','std');

fh = figure(1); clf
plot(Torque.times,Ang,'LineWidth',2)
title('Click to continue')
hold on
for ii = 1:length(pts)
    xline(Torque.times(pts(ii)),'LineWidth',2,'Color','r')
end
ylabel('Angle (°)')

yyaxis right
plot(Torque.times,Tor,'LineWidth',2)
title('Click to continue')
ylabel('Tor (Nm)')
xlabel('Time (s)')
waitforbuttonpress

%% Poly for passive torque
polyp = Auto_RestTorque(Tor,Ang,Frequency);

%% Torque values according to the phases
% c = 0.5 * Frequency; %just to get rid of possible peaks due to people not relaxing before the end, only half second
% Rest_Sho = round( mean(Tor(1: pts(1) -c ) ) ,2 );
% Rest_Long = round( mean( Tor(pts(2) : pts(3) ) ) ,2);

%remove outliers as maybe people did not rest completely/artifacts
Rest_Sho = mean(rmoutliers(round( Tor(1: pts(1) ), 2 ) ));
Rest_Long = mean(rmoutliers(round( Tor(pts(2) : pts(3) ), 2 ) ));

%% Rotation stuff
Tor_rot = round( Tor(pts(1):pts(2)), 2);
Ang_rot = round( Ang(pts(1):pts(2)), 2);

%polynomial to get rid of rest torque and just scale active torque
[Pas_Tor_Rot,Tor_rot] = Normalize_Torque(Tor_rot,Ang_rot,polyp);

%% Ask scaling and duration of rotation
% it's better to scale here as I can scale JUST the dynamic/movement part
% and not the resting torque
ramp = (pts(2)-pts(1)) /Frequency;

% clf
% prompt = {'Enter Ramp duration(s)','Enter Scaling Factor (0-1)'};
% dlgtitle = 'Ramp design';
% dims = [1 35];
% definput = {num2str(ramp),'0.75'};
%
% answer = inputdlg(prompt,dlgtitle,dims,definput);

% ramp = str2double(answer{1});
% factorScaling = str2double(answer{2});
[ramp, factorScaling, ContractionType] = dlg_SetTrace(ramp);

%% Create rotation phase according to the input(s)
nPoint_contraction = Frequency * ramp;
Tor_rot = Tor_rot .* factorScaling; %scale net torque during rotation
Tor_rot =  Tor_rot + Pas_Tor_Rot; %re-add passive to scaled active torque

%interpolate to trace points according to the frequency and for having a
%good "shape" in the xy view
Tor_Con = interp1( linspace(0,1,length(Tor_rot)), Tor_rot, linspace(0,1,nPoint_contraction));
Tor_Ecc = flip(Tor_Con);

%% Assemble the file

Time = [-1 1.99 linspace(2,2+ramp,nPoint_contraction) 2+ramp 8 9 10]; %time is equal for both
for ii =  1 : size(Tor_rot,2)
    Concentric(:,ii) = [Rest_Sho Rest_Sho Tor_Con(:,ii)' Rest_Long Rest_Long Rest_Sho Rest_Sho]';
    Eccentric(:,ii) = [Rest_Long Rest_Long Tor_Ecc(:,ii)' Rest_Sho Rest_Sho Rest_Long Rest_Long]';
end

% clf
% plot(Time,Concentric,'LineWidth',2);
% hold on
% plot(Time,Eccentric,'LineWidth',2);
% legend('Concentric','Eccentric')

%% Save Matrix txts for Spike2 plotting
%
if size(factorScaling,2) > 3
    WeekTxt = '_P1_';
else WeekTxt = '_P2_';
end

Matrix(:,1) = Time;
clf

if strcmp(ContractionType,'Concentric')

    %[fname,pth] = uiputfile(name,'Where do you wanna save the CONCENTRIC Trace?'); % Type in name of the file.
    pth = uigetdir(pth,'Where do you wanna save the CONCENTRIC Traces?');

    %save each intensity trace in a separated file
    for ii = 1 : size(Concentric,2)
        Matrix(:,2) = Concentric(:,ii);

        fname = ['/ConMatch' WeekTxt num2str(factorScaling(ii)*100) '.txt'];
        if pth ~= 0
            writematrix(Matrix,[pth fname]);
        end

        plot(Time,Concentric(:,ii),'LineWidth',2);
        hold on
    end
else
    %get the directory where to save the traces
    pth = uigetdir(pth,'Where do you wanna save the ECCENTRIC Traces?');

    for ii = 1 : size(Eccentric,2)

        Matrix(:,2) = Eccentric(:,ii);

        fname = ['/EccMatch' WeekTxt num2str(factorScaling(ii)*100) '.txt'];

        if pth ~= 0
            writematrix(Matrix,[pth fname]);
        end
        plot(Time,Eccentric(:,ii),'LineWidth',2);
        hold on
    end

end
