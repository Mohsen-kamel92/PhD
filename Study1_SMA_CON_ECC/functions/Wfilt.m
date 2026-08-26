function yF = Wfilt(y,fc,type,fs) % Take signal y, filter it with cutoff fc, 
% using filter type'string input', assuming sampling rate fs, and return filtered signal yF


%*Filter signal with a 2nd-order Butterworth filter*
%Wfilt(signal, cut-off_frequency_in_Hz, filter_type, sampling_frequency)
%*Requires Wcu function.

% Author:
% BJ Raiteri, April 2021, if you find errors pls email brent.raiteri@rub.de
% tested in R2019b

% Change signal to column vector
if isrow(y)
    y = y';
end

% Remove NaNs from signal
if sum(isnan(y)) > 0 
    tf = isnan(y); % detect NaNs
    idx = 1:numel(y); % creat an index of y
    y(tf) = interp1(idx(~tf),y(~tf),idx(tf)); %Replace NaNs using linear interpolation between valid neighbors.
    if ~isempty(isnan(y)) % If NaNs exist at beginning or end, interpolation cannot estimate them, so cut the tiles outside the interpolation window
        cut = find(isnan(y),1,'first');
        y = y(1:cut-1); % values are save of NaNs now, because if NaNs exist, it will propagate and filter will result in NaNs in the end
    end      
end
% low pass: keep frequencies & ermove high frequencies
% Change filter type based on type

if strcmp(type,'low')
    [cu,order] = Wcu(4, fc, 'low', fs);   % 1st input = 4 for 2nd-order filter
    [B,A] = butter(order,cu/(fs/2),'low'); % cutoff here is a fraction of Nyquist, so
    % (The highest frequency you can correctly represent is half the
    % sampling rate, So you get exactly two samples per cycle in case
    % reassigned 2k to 1k).
    % For coefficients [B,A], based on a 2nd-order Butterworth filter at 20
    % Hz, butterworth formula computes the exact coefficients that satisfy
    % the filter order inputs in time-steps
elseif strcmp(type,'high')
    [cu,order] = Wcu(4, fc, 'high', fs);   % 1st input = 4 for 2nd-order filter
    [B,A] = butter(order,cu/(fs/2),'high');
elseif strcmp(type,'bandpass')
    if size(fc,2) ~= 2
        error('Provide two cut-off frequencies for a bandpass filter.')
    end
    cu = fc;
    order = 2;  
    [B,A] = butter(order,cu/(fs/2),'bandpass');
else
    error('Enter a valid filter type [low high bandpass stop].')
end

% Filter signal 
yF = filtfilt(B,A,y);