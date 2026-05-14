function [ varargout ] = Normalize_Torque (TorSig,AngSig,p)

%% From passive data fitted, net torque at that specific angle is calculated and subtracted
%  Author: Paolo Tecchio (Paolo.Tecchio@rub.de)


%% check arguments
minArgs=1;  
maxArgs=2;
nargoutchk(minArgs,maxArgs);

%disp("You requested " + nargout + " outputs.")
varargout = cell(nargout,1);

%initialize variables
Pas_Tor = zeros(length(AngSig),1);
Net_Tor = zeros(length(AngSig),1);

%Passive torque & Net torque calculation
Pas_Tor = polyval(p,AngSig);
Net_Tor = TorSig-Pas_Tor;

%% Outputs
%first output is the passive, second is the Net 
varargout{1} = Pas_Tor;
varargout{2} = Net_Tor;
end