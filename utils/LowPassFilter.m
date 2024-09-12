%% LowPassFilter
% 
% Apply a low-pass filter and downsample the data.
% 
% Description:
%   This function applies a low-pass Butterworth filter and then downsamples the input data. 
%   The filter is designed based on the specified cutoff frequency (fc( and sampling frequency (fs). 
% 
% Inputs:
%   X      - [channels, time points] matrix of data to be filtered and downsampled.
%   ds     - Downsampling factor (integer). The output data will contain every `ds`-th sample.
%   fs     - Sampling frequency of the data in Hz (optional, default is 1000 Hz).
%   fc     - Cutoff frequency of the low-pass filter in Hz (optional, default is 100 Hz).
% 
% Outputs:
%   Y      - [channels, downsampled time points] matrix of filtered and downsampled data.
% 
% Alberto Liardi, 2024

function [Y] = LowPassFilter(X,ds,fs,fc)

    % default parameters
    if ~exist('fs', 'var'), fs = 1000; end % sampling frequency 
    if ~exist('fc', 'var'), fc = 100; end % frequency cut off
    
    % getting coefficient of the filter
    [b, a] = butter(2, fc/(fs/2));
    
    % filtering along rows
    Y = filter(b,a,X,[],2); 
    
    % downsampling
    Y = Y(1:ds:end);

end