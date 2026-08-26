function signalVelocity = Wvel(y,x)
%*Calculate velocity based on 2x*
%Wvel(signal, time)

% Author:
% BJ Raiteri, April 2021, if you find errors pls email brent.raiteri@rub.de
% tested in R2022a

if iscolumn(y)
    y = y';
end

% If x is a scalar then create a time vector
if isscalar(x)
    x = 1/x:1/x:1/x*length(y);
elseif iscolumn(x)
    x = x';
end

signalVelocity = zeros(size(y));

for kk = 2:length(y)-1

    signalVelocity(1,kk) = (y(1,kk+1)-y(1,kk-1))/(x(1,kk+1)-x(1,kk-1));

end

end