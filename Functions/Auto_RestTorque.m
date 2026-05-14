function  Passive_p  = Auto_RestTorque (Tor, Ang, Frequency)
%%Automatic detection of rest torque at two positions
%% Filtering w
w= abs(diff(Ang) ./ (1/Frequency));
%F_Spectrum(w)

CutOff = 4;
Filt_order = 6;

%correct order and cut off
[cu,order] = Wcu(Filt_order,CutOff,'low',Frequency);
[b,a] = butter(order, cu/(0.5*Frequency),'low');

w = filtfilt(b,a,w);
%% Find points where the velocity is close to the minimum (zero) with a threshold error
threshold = 2;
locs = find(w <= (abs(min(w)) + threshold));


% split value according to the two positions
Tor_L = nonzeros(Tor(locs) .* (Tor(locs) < mean(Tor(locs))));
Ang_L = nonzeros(Ang(locs) .* (Tor(locs) < mean(Tor(locs))));

Tor_S = nonzeros(Tor(locs) .* (Tor(locs) > mean(Tor(locs))));
Ang_S = nonzeros(Ang(locs) .* (Tor(locs) > mean(Tor(locs))));

% clf
% plot(Ang_L,Tor_L,'r*')
% hold on
% plot(Ang_S,Tor_S,'g*')

%% Clear outliers as potential value could be close to the very end of the rotations 
% the subject did not relax / drop the foot at end
%Outliers are defined as elements more than three scaled MAD from the median. 
% The scaled MAD is defined as c*median(abs(A-median(A))), where c=-1/(sqrt(2)*erfcinv(3/2)).
[Tor_L, idx] = rmoutliers(Tor_L,'median');
Ang_L =  nonzeros(Ang_L .* ~idx)';

[Tor_S, idx] = rmoutliers(Tor_S,'median');
Ang_S =  nonzeros(Ang_S .* ~idx)';

% plot(Ang_L,Tor_L,'g*')
% hold on
% plot(Ang_S,Tor_S,'g*')

%% Passive tor during rotation back

%threshold = 2;
w = abs(abs(w)-40); %get passive torque during returning rotation
locs = find(w <= threshold);
%% Poly fit and save it
%fit 1st poly
Passive_p = polyfit([Ang_S Ang_L Ang(locs)'],[Tor_S' Tor_L' Tor(locs)'],1);
% 
% Ang = round(min(Ang)-2 : 1 : max(Ang)+2);
% Tor = polyval(Passive_p,Ang);
% 
% title('Rest fit torque')
% hold on
% plot(Ang,Tor,'--','LineWidth',4,'color',[0 0.21 0.38]);
% ylabel('Torque (Nm)')
% xlabel('Crank Foot Angle (Degree)');
% set(gca,'FontSize',16)
% box off
