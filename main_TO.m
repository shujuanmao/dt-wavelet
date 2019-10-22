%%    An example of calculating dt using My_Wxspectrum.m
%     modify ParameterChoice in line 44 to test with different sets of parameters        
%         
%     Created: Jul., 2019, 
%     Shujuan Mao (maos@mit.edu) and Aurélien Mordret (mordret@mit.edu)
%
%%
%%    Reference: 
%           S.Mao, A.Mordret, M.Campillo, H.Fang, R.D. van der Hilst,(2019),
%           On the Measurement of Seismic Travel-Time Changes in the
%           Time-Frequency Domain with Wavelet Cross-Spectrum Analysis,
%           GJI, In Review.
%
%%
%%    Copyright (c) 2019, Shujuan Mao and Aurélien Mordret, covered by MIT License.
% 
%     Permission is hereby granted, free of charge, to any person obtaining a copy
%     of this software and associated documentation files (the "Software"), to deal
%     in the Software without restriction, including without limitation the rights
%     to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%     copies of the Software, and to permit persons to whom the Software is
%     furnished to do so, subject to the following conditions:
% 
%     The above copyright notice and this permission notice shall be included in all
%     copies or substantial portions of the Software.
% 
%     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%     IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%     AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%     LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%     SOFTWARE.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

close all;
clear;
%% Read in synthetic data
% Difference between current and reference: 0.05% dv/v homogeneously
% (for more details see Section 3.1 in the reference.)
load('synthetic_dvov_0.05percent.mat');


%% Parameters for My_Wxspectrum.m
%  See notes in My_Wxspectrum.m for details of each parameter

ParameterChoice = 1;  % 1 or 2, to try with two sets of parameters
                      
if (ParameterChoice == 1)
    wname = 'amor';
    WaveletParameters = [];
    FrequencyLimits = [0.5,5];  % in Hz
    SmoothingFlag = 1;
    NumScalesToSmooth = 3;
    DegTimeToSmooth = 0.25;
    VoicesPerOctave = 10;
    ExtendSigFlag = 1;
    
elseif (ParameterChoice == 2)
    wname = 'morse';
    WaveletParameters = [3,80];
    FrequencyLimits = [0.2,8];  % in Hz
    SmoothingFlag = 0;
    NumScalesToSmooth = [];
    DegTimeToSmooth = [];
    VoicesPerOctave = 16;
    ExtendSigFlag = 1;
    
else
    fprintf('Wrong value for ParameterChoice!\n')
end


%% Calculating the time-frequency distribution of dt
[WXspec,WXdt,WXamp,Wcoh,Freq,Coi] = My_Wxspectrum_TO(ori_waveform,new_waveform,Fs,wname,...
    WaveletParameters,FrequencyLimits,SmoothingFlag,NumScalesToSmooth,DegTimeToSmooth,...
    VoicesPerOctave, ExtendSigFlag);

% If not need to calculate Wcoherence:
% [WXspec,WXdt,WXamp,Freq,Coi] = My_Wxspectrum_TO(ori_waveform,new_waveform,Fs,wname,...
%     WaveletParameters,FrequencyLimits,SmoothingFlag,NumScalesToSmooth,DegTimeToSmooth,...
%     VoicesPerOctave, ExtendSigFlag);


%% Plotting
figure()
subplot(221)
hp = pcolor(time,Freq, WXdt);
set(hp,'edgecolor','none');
colormap jet; 
colorbar;
caxis(0.05.*[-1,1]);
hold on;
ax = gca;
plot(ax, time, Coi,'w--','linewidth',2);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Time difference');
set(gca,'fontsize',13);

subplot(222)
hp = pcolor(time,Freq, angle(WXspec));
set(hp,'edgecolor','none');
colormap jet; 
colorbar;
caxis(1.*[-1,1]);
hold on;
ax = gca;
plot(ax, time, Coi,'w--','linewidth',2);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Phase difference');
set(gca,'fontsize',13);

subplot(223)
hp = pcolor(time,Freq, log(WXamp));
set(hp,'edgecolor','none');
colormap jet; 
colorbar;
hold on;
ax = gca;
plot(ax, time, Coi,'w--','linewidth',2);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('log of local power');
set(gca,'fontsize',13);

subplot(224)
hp = pcolor(time,Freq, Wcoh);
set(hp,'edgecolor','none');
colormap jet; 
colorbar;
caxis([0.9,1]);
hold on;
ax = gca;
plot(ax, time, Coi,'w--','linewidth',2);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Wavelet coherence');
set(gca,'fontsize',13);

