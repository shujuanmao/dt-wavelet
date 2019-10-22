function varargout= My_Wxspectrum_TO(x_reference,x_current,fs,wname,WaveletParameters,...
    FrequencyLimits,SmoothingFlag,NumScalesToSmooth,DegTimeToSmooth,...
    VoicesPerOctave, ExtendSigFlag)
%% My_Wxspectrum_TS: 
%   Using wavelet cross-spectrum to calculate the lapse-time- and frequency-dependent 
% travel-time changes between two time series. !!!! Matlab R2018 (or 2019) and 
% the WAVELET Toolbox are required for this code. !!!!
%
%
%% USAGE
%  [WXspec,WXdt,WXamp,Wcoh,Freq,Coi]= My_Wxspectrum_TO(varargin)
%  OR:
%  [WXspec,WXdt,WXamp,Freq,Coi]= My_Wxspectrum_TO(varargin)
%  where varargin = [x,y,fs,wname,WaveletParameters,FrequencyLimits,...
%                    SmoothingFlag,NumScalesToSmooth,DegTimeToSmooth,...
%                    VoicesPerOctave,ExtendSigFlag]
%     
%  
%  Note: Usage2 does not provide "Wcoher" in the output, and the calulation 
%   is faster compared with the Usage1 when the SmoothingFlag is truned off.
%
%% Input
%    x_reference,x_current: Two vectors, reference and current time series.
%    fs: Positive scalar, sampling frequency
%    wname: The type of the wavelet. Avalible options include 
%           1) 'amor' for Morlet wavelet; 2) 'morse' for Morse wavelet; 
%           and 3) 'bump' for bump wavelet. Recommanded: wname = 'amor'.
%    WaveletParameters: For wname = 'morse', set WaveletParameters = [gamma, P2] 
%           or as [] to use default value [3,60]; for the other two types of
%           wavelets, set WaveletParameters = [].
%    FrequencyLimits: Range of frequency to calculate CWT (a two-element vector)
%    SmoothingFlag: 0 or 1, corresponding to whether or not to smooth CWT results.
%                   Recommanded: SmoothingFlag = 0 .
%                   If specified as 1, a boxcar window will be applied in
%                   the scale direction and a Gaussian window in the time
%                   direction for smoothing. !!!! This type of smoothing is only 
%                   appropriate for Morlet wavelet (wname = 'amor'). !!!!
%    NumScalesToSmooth: Positive integer, indicting the length of boxcar window. 
%    DegTimeToSmooth: Positive scalar,indicating the length of the Gaussina window. 
%                     If set DegTimeToSmooth = [], default value of 0.25 will be used; 
%                     larger values mean more smoothing, smaller mean less smoothing.
%    (NumScalesToSmooth = 1 and DegTimeToSmooth = 0 means no smoothing in
%    either directions.)
%    VoicesPerOctave: Even integer from 4 to 48, indicates how fine the frequency 
%                     is discretized. Recommanded to be no less than 10.
%    ExtendSigFlag: 0 or 1, corresponding to whether or not to extend the signal
%                   symmetrically to mitigate boundary effects.
%
%
%% OUTPUT
%    WXspec: complex-valued matrix, the wavelet cross-spectrum
%    WXdt: matrix of time difference and phase difference, respectively
%                  between the two input time series in time-frequency domain.
%                  !!!! This WXdt is obtained with wrapped phase difference.
%                  If needed, the user can also produce time difference with 
%                  unwrapped phase from angle(WXspec). !!!!
%    WXamp: matrix of amplitude product of two CWT in time-frequency domain
%    Wcoher: matrix of wavelet coherence 
%    Freq: vector of frequenies used in CWT, in Hz
%    Coi: Cone of incluence, indicating areas affected by edge effects.
%
%
%% EXAMPLE:
%    [WXspec,WXdt,WXamp,Wcoher,Freq,Coi]= My_Wxspectrum(current,reference, ...
%                                          50,'amor',[],[0.5,5],1,3,0.25,10,1)
%  OR
%    [WXspec,WXdt,WXamp,Freq,Coi]= My_Wxspectrum(current,reference, ...
%                                          50,'amor',[],[0.5,5],1,3,0.25,10,1)
%
%%    Authors: Shujuan Mao (maos@mit.edu) and Aurélien Mordret (mordret@mit.edu)
%     Created: Aug., 2018
%     Updated: Jul., 2019
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

    % check the number of outputs
    if (nargout ~=5) && (nargout ~=6)
        error('Wrong number of outputs. Output should be formatted in either [WXspec,WXdt,WXamp,Wcoher,Freq,Coi] or [WXspec,WXdt,WXamp,Freq,Coi].');
    end
    
    
    %Check the inputs
    if all(size(x_current) ~= 1) || all(size(x_reference) ~= 1)
	    error('x and y must be vectors, not matrix.')
    end
    nx = numel(x_current); 
    ny = numel(x_reference);
    if (~isequal(nx,ny) || numel(x_current) < 4)
        error(message('Wavelet:FunctionInput:EqualLengthInput'));
    end 
    
    if (VoicesPerOctave > 48) || (VoicesPerOctave < 4) || (mod(VoicesPerOctave,2) ~= 0)
        error('VoicesPerOctave must be an even integer between 4 and 48.')
    end
    

    % Form signals as row vectors
    x_reference = x_reference(:)';  %%% reference
    x_current = x_current(:)';  %%% current
    

    dt = 1/fs;

    nv = VoicesPerOctave;

    if (isempty(DegTimeToSmooth))
        nt = 0.25;
    else
        nt = DegTimeToSmooth;
    end

    if (isempty(NumScalesToSmooth))
        ns = 3;
    else
        ns = NumScalesToSmooth;
    end


    if (strcmp(wname, 'morse'))
        if (isempty(WaveletParameters))
            WaveletParameters = [3,60];   
        end
        
        [cwt_reference,Freq,Coi,fb] = cwt(x_reference,wname,fs,'ExtendSignal',ExtendSigFlag,'VoicesPerOctave',nv,...
            'WaveletParameters',WaveletParameters, 'FrequencyLimits', FrequencyLimits);
        cwt_current = cwt(x_current,wname,fs,'ExtendSignal',ExtendSigFlag,'VoicesPerOctave',nv,...
            'WaveletParameters',WaveletParameters, 'FrequencyLimits', FrequencyLimits);
        
       
    elseif (strcmp(wname, 'amor') || strcmp(wname, 'bump'))
        [cwt_reference,Freq,Coi,fb] = cwt(x_reference,wname,fs,'ExtendSignal',ExtendSigFlag,'VoicesPerOctave',nv, ...
            'FrequencyLimits', FrequencyLimits);
        cwt_current = cwt(x_current,wname,fs,'ExtendSignal',ExtendSigFlag,'VoicesPerOctave',nv, ...
            'FrequencyLimits', FrequencyLimits);
        
      
    else
        error('Unvalid input of the wavelet name! Avalible wavelets: morse, amor, bump');   
    end
    
    scales = fb.scales;
    scales = scales';
    invscales = 1./scales;
    invscales = repmat(invscales,1,nx);

    if ((~SmoothingFlag)||(ns == 1 && nt == 0))
        fprintf('Without Smoothing\n')
        %%% Without smoothing
        crossCFS = cwt_reference.*conj(cwt_current);
        WXamp = abs(crossCFS);
        WXspec = crossCFS;
        
        %%%% For calculating coherence
        if (nargout == 6)
            cfs1 = smoothCFS(invscales.*abs(cwt_current).^2,scales,dt,ns,nt);
            cfs2 = smoothCFS(invscales.*abs(cwt_reference).^2,scales,dt,ns,nt);
            crossCFS = smoothCFS(invscales.*crossCFS,scales,dt,ns,nt);
            Wcoh = abs(crossCFS).^2./(cfs1.*cfs2);
        end


    else
        fprintf('With Smoothing\n')
        % With smoothing
        cfs1 = smoothCFS(invscales.*abs(cwt_current).^2,scales,dt,ns,nt);
        cfs2 = smoothCFS(invscales.*abs(cwt_reference).^2,scales,dt,ns,nt);
        crossCFS = cwt_reference.*conj(cwt_current);
        WXamp = abs(crossCFS);

        crossCFS = smoothCFS(invscales.*crossCFS,scales,dt,ns,nt);
        WXspec = crossCFS./(sqrt(cfs1).*sqrt(cfs2));
        Wcoh = abs(crossCFS).^2./(cfs1.*cfs2);
    end

    WXangle = angle(WXspec);
    WXdt = WXangle./repmat(2.*pi.*Freq,1,nx);
    
   
    varargout{1} = WXspec;
    varargout{2} = WXdt;
    varargout{3} = WXamp;
    if (nargout == 5)
        varargout{4} = Freq;
        varargout{5} = Coi;
    else
        varargout{4} = Wcoh;
        varargout{5} = Freq;
        varargout{6} = Coi;
    end
    
end

%% A particular choice of smoothing
function cfs = smoothCFS(cfs,scales,dt,ns,nt)

    N = size(cfs,2);
    npad = 2.^nextpow2(N);
    omega = 1:fix(npad/2);
    omega = omega.*((2*pi)/npad);
    omega = [0., omega, -omega(fix((npad-1)/2):-1:1)];

    % Normalize scales by DT because we are not including DT in the
    % angular frequencies here. The smoothing is done by multiplication in
    % the Fourier domain
    normscales = scales./dt;
    for kk = 1:size(cfs,1)
        F = exp(-nt*(normscales(kk)^2)*omega.^2);
        smooth = ifft(F.*fft(cfs(kk,:),npad));

        cfs(kk,:)=smooth(1:N);
    end
    
    % Convolve the coefficients with a moving average smoothing filter across scales
    H = 1/ns*ones(ns,1);
    cfs = conv2(cfs,H,'same');

end



