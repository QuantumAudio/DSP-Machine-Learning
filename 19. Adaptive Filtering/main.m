%% 1.1: Adaptive cancellation of power-frequency interference
s = load('s60.dat'); % true ECG
x = load('x60.dat'); % ECG + 60 Hz + 180 Hz interference
y = load('y60.dat'); % 60 Hz + 180 Hz reference

fs = 1000; % sampling rate;

%% 1.1.a: Adaptive Wiener Filter
N = length(x);
M = 50;
mu = 1e-5;
h = zeros(M+1,1);   % initialize empty filter values
w = zeros(M+1,1);   % initialize internal state buffer
% LMS algorithm
for n = 1:N
    w(1) = y(n);    % load internal state buffer
    xHat(n) = h'*w; % use h and y(n) to estimate of x(n)
    e(n) = x(n) - xHat(n);  % error output
    h = h + 2*mu*e(n)*w;    % update filter weights
    w = [w(1); w(1:end-1)]; % update buffer
end

%
t = 0:1/fs:(2 - 1/fs);       % time index

figure
plot(t,x), title('x = ECG + 60&180 Hz interference'),xlabel('t (sec)')
axis([0, 2, -6, 6]);
%
figure
plot(t,y), title('y = 60&180 Hz reference'),xlabel('t (sec)')
axis([0, 2, -6, 6]);
%
figure
plot(t,e,t,s,'r'), title('e = estimated ECG, \mu = 1.0e-5'),xlabel('t (sec)')
axis([0, 2, -6, 6]), legend('LMS algorithm','true ECG','location','southeast')
%
%% 1.1.b: Repeat 1.1.a with mu = 2e-5
mu = 2e-5;

h = zeros(M+1,1);   % initialize empty filter values
w = zeros(M+1,1);   % initialize internal state buffer
% LMS algorithm
for n = 1:N
    w(1) = y(n);    % load internal state buffer
    xHat(n) = h'*w; % use h and y(n) to estimate of x(n)
    e(n) = x(n) - xHat(n);  % error output
    h = h + 2*mu*e(n)*w;    % update filter weights
    w = [w(1); w(1:end-1)]; % update buffer
end


% Note the slowing down of the LMS learning speed
figure
plot(t,e,t,s,'r'), title('e = estimated ECG, \mu = 2.0e-5'),xlabel('t (sec)')
axis([0, 2, -6, 6]), legend('LMS algorithm','true ECG','location','southeast')

%
% Experiment with different mus and determine the largets mu we can use
% without the LMS algorithm diverging
% mu = 1.0e-3; % signal is not exact but still useful
mu = 3.0e-3; % just barely converges
% mu = 3.1e-3; % mu diverges

h = zeros(M+1,1);   % initialize empty filter values
w = zeros(M+1,1);   % initialize internal state buffer
% LMS algorithm
for n = 1:N
    w(1) = y(n);    % load internal state buffer
    xHat(n) = h'*w; % use h and y(n) to estimate of x(n)
    e(n) = x(n) - xHat(n);  % error output
    h = h + 2*mu*e(n)*w;    % update filter weights
    w = [w(1); w(1:end-1)]; % update buffer
end

figure
plot(t,e,t,s,'r'), title('e = estimated ECG, \mu = 3.0e-3'),xlabel('t (sec)'),
axis([0, 2, -6, 6]), legend('LMS algorithm','true ECG','location','southeast')

%
%% 1.1.c: Repeat part 1.1.a with RLS algorithm, la = 0.99, delta = 0.1
la = 0.99; delta = 0.1; P = speye(M+1)/delta;
w = zeros(M+1,1); h = zeros(M+1,1); 

clear k; clear k1; clear e; clear eEst; clear xEst; 

% RLS algorithm
for n = 1:N
    w(1) = y(n);
    k = (1/la)*P*w;
    v = k'*w;
    mu = 1/(1 + v);
    k1 = mu*k;
    P = (1/la)*P - k1*k';
    P = (1/2)*(P + P');
    xEst = h'*w;
    eEst = x(n) - xEst;
    e(n) = mu*eEst;
    xEst = x(n) - e(n);
    h = h + eEst*k1;
    w = [w(1); w(1:end-1)];
end

figure
plot(t,e,t,s,'r'), title('e = estimated ECG, \lambda = 0.990, \delta = 0.1'),xlabel('t (sec)'),
axis([0, 2, -6, 6]), legend('RLS algorithm','true ECG','location','southeast')
%
% Repeat with delta = 10
la = 0.99; delta = 10; P = speye(M+1)/delta;
w = zeros(M+1,1); h = zeros(M+1,1); 

clear k; clear k1; clear e; clear eEst; clear xEst; 

% RLS algorithm
for n = 1:N
    w(1) = y(n);
    k = (1/la)*P*w;
    v = k'*w;
    mu = 1/(1 + v);
    k1 = mu*k;
    P = (1/la)*P - k1*k';
    P = (1/2)*(P + P');
    xEst = h'*w;
    eEst = x(n) - xEst;
    e(n) = mu*eEst;
    xEst = x(n) - e(n);
    h = h + eEst*k1;
    w = [w(1); w(1:end-1)];
end

figure
plot(t,e,t,s,'r'), title('e = estimated ECG, \lambda = 0.990, \delta = 10'),xlabel('t (sec)'),
axis([0, 2, -6, 6]), legend('RLS algorithm','true ECG','location','southeast')
%
% Repeat with delta = 100

la = 0.99; delta = 100; P = speye(M+1)/delta;
w = zeros(M+1,1); h = zeros(M+1,1); 

clear k; clear k1; clear e; clear eEst; clear xEst; 

% RLS algorithm
for n = 1:N
    w(1) = y(n);
    k = (1/la)*P*w;
    v = k'*w;
    mu = 1/(1 + v);
    k1 = mu*k;
    P = (1/la)*P - k1*k';
    P = (1/2)*(P + P');
    xEst = h'*w;
    eEst = x(n) - xEst;
    e(n) = mu*eEst;
    xEst = x(n) - e(n);
    h = h + eEst*k1;
    w = [w(1); w(1:end-1)];
end

figure
plot(t,e,t,s,'r'), title('e = estimated ECG, \lambda = 0.990, \delta = 100'),xlabel('t (sec)'),
axis([0, 2, -6, 6]), legend('RLS algorithm','true ECG','location','southeast')

%
% Repeat with delta = 1000
la = 0.99; delta = 1000; P = speye(M+1)/delta;
w = zeros(M+1,1); h = zeros(M+1,1); 

clear k; clear k1; clear e; clear eEst; clear xEst; 

% RLS algorithm
for n = 1:N
    w(1) = y(n);
    k = (1/la)*P*w;
    v = k'*w;
    mu = 1/(1 + v);
    k1 = mu*k;
    P = (1/la)*P - k1*k';
    P = (1/2)*(P + P');
    xEst = h'*w;
    eEst = x(n) - xEst;
    e(n) = mu*eEst;
    xEst = x(n) - e(n);
    h = h + eEst*k1;
    w = [w(1); w(1:end-1)];
end

figure
plot(t,e,t,s,'r'), title('e = estimated ECG, \lambda = 0.990, \delta = 1000'),xlabel('t (sec)'),
axis([0, 2, -6, 6]), legend('RLS algorithm','true ECG','location','southeast')

%% 1.1.d: Using a double notch filter to remove 60 & 180 Hz interference
[b1, a1, beta1] = parmeq(1, 0, 1/sqrt(2), 2*pi*60/fs, 2*pi*3/fs);
[b2, a2, beta2] = parmeq(1, 0, 1/sqrt(2), 2*pi*180/fs, 2*pi*3/fs);

y1 = filter(b1, a1, x);
y2 = filter(b2, a2, y1);


figure
plot(t,y1,t,s,'r'), title('filtered ECG - single notch'),xlabel('t (sec)'),
axis([0, 2, -6, 6]), legend('filtered','true ECG','location','southeast')

figure
plot(t,y2,t,s,'r'), title('filtered ECG - double notch'),xlabel('t (sec)'),
axis([0, 2, -6, 6]), legend('filtered','true ECG','location','southeast')

%% 1.1.e: Periodogram spectra
% Spectrum of x
f0 = 60; % Dominant frequency in [Hz]

[Px,f] = pwelch(x,[],[],[],fs);     % power spectrum estimate
Pdb = 10*log10(abs(Px));            % convert to dB

figure
plot(f/f0,Pdb); xlabel('f/f_0'), ylabel('dB'), 
axis([0, 8, -100, 5]),
grid on, title('|S_{xx}(f)|, ECG + 60&180 Hz noise')

% Spectrum of y
[Py,f] = pwelch(y,[],[],[],fs);     % power spectrum estimate
Pdb = 10*log10(abs(Py));            % convert to dB

figure
plot(f/f0,Pdb); xlabel('f/f_0'), ylabel('dB'), 
axis([0, 8, -100, 5]),
grid on, title('|S_{yy}(f)|, noise reference y')

%
% Spectrum of e

% Repeat of 2.1.a (e had been overwritten)
N = length(x);
M = 50;
mu = 1e-4;
h = zeros(M+1,1);   % initialize empty filter values
w = zeros(M+1,1);   % initialize internal state buffer
% LMS algorithm
for n = 1:N
    w(1) = y(n);    % load internal state buffer
    xHat(n) = h'*w; % use h and y(n) to estimate of x(n)
    e(n) = x(n) - xHat(n);  % error output
    h = h + 2*mu*e(n)*w;    % update filter weights
    w = [w(1); w(1:end-1)]; % update buffer
end

[Pe,f] = pwelch(e,[],[],[],fs);     % power spectrum estimate
Pdb = 10*log10(abs(Pe));            % convert to dB

figure
plot(f/f0,Pdb); xlabel('f/f_0'), ylabel('dB'), 
axis([0, 8, -100, 5]),
grid on, title('|S_{ee}(f)|, LMS-recovered ECG')

%
% Spectrum of double-notch filtered signal y2
[Py2,f] = pwelch(y2,[],[],[],fs);     % power spectrum estimate
Pdb = 10*log10(abs(Py2));            % convert to dB

figure
plot(f/f0,Pdb); xlabel('f/f_0'), ylabel('dB'), 
axis([0, 8, -100, 5]),
grid on, title('|S_{nf}(f)|, double-notch filtered')

%% 1.2: Adaptive separation of maternal and fetal ECGs

x = zmean(load('M1+B.dat'));    % maternal + fetal ECG
y = zmean(load('M2.dat'));      % maternal reference ECG from thorax

% fs = 1000; (same as previously defined)

[x, xMean] = zmean(x'); % Make x zero-mean and the mean value
[y, yMean] = zmean(y'); % Make y zero-mean and the mean value 

t = (1:length(x))/fs;       % time index scaled for plotting

figure
plot(t,x), title('mother + baby'),xlabel('t (sec)')
axis([0, 9.5, -1.1, 3.1])

figure
plot(t,y,'m'), title('mother reference'),xlabel('t (sec)')
axis([0, 9.5, -1.1, 3.1])

% Run LMS algorithm on x and y
N = length(x);
M = 90;              % Note: M < 50 did not extract enough interference
% and M > 150 started to smooth desired transients 
mu = 1e-3;
h = zeros(M+1,1);   % initialize empty filter values
w = zeros(M+1,1);   % initialize internal state buffer
% LMS algorithm
for n = 1:N
    w(1) = y(n);    % load internal state buffer
    xHat(n) = h'*w; % use h and y(n) to estimate of x(n)
    e(n) = x(n) - xHat(n);  % error output
    h = h + 2*mu*e(n)*w;    % update filter weights
    w = [w(1); w(1:end-1)]; % update buffer
end

figure
plot(t,e,'r'), title('baby ECG, LMS \mu = 0.0010'),xlabel('t (sec)')
axis([0, 9.5, -1, 1]);


% Repeat with mu = 1e-4

% Run LMS algorithm on x and y
N = length(x);
M = 90;              % Note: M = 1 even looked fine
mu = 1e-4;
h = zeros(M+1,1);   % initialize empty filter values
w = zeros(M+1,1);   % initialize internal state buffer
% LMS algorithm
for n = 1:N
    w(1) = y(n);    % load internal state buffer
    xHat(n) = h'*w; % use h and y(n) to estimate of x(n)
    e(n) = x(n) - xHat(n);  % error output
    h = h + 2*mu*e(n)*w;    % update filter weights
    w = [w(1); w(1:end-1)]; % update buffer
end

t = (1:length(x))/fs;       % time index scaled for plotting

figure
plot(t,e,'r'), title('baby ECG, LMS \mu = 0.0001'),xlabel('t (sec)')
axis([0, 9.5, -1, 1]);
 
%% 1.3: Filtering a noisy fetal heartbeat sound

[x,fs] = audioread('B+noise.wav');      % noisy heartbeat sound
[y,fs] = audioread('B-filtered.wav');   % filtered signal

% Note: fs = 22.05 [kHz]

t = (1:length(x))/fs;       % time index scaled for plotting

figure
plot(t,x), title('noisy heartbeat'),xlabel('t (sec)')
axis([0, 6.5, -0.5, 1])

% 1.3.a
[x, xMean] = zmean(x'); % Make x zero-mean and the mean value
[y, yMean] = zmean(y'); % Make y zero-mean and the mean value 


% Spectrum of x
[Px,f] = pwelch(x,[],[],[],fs);     % power spectrum estimate
Pdb = 10*log10(abs(Px));            % convert to dB


figure          % Note: There seems to be an offset of about -45 dB
plot(f/1000,Pdb); xlabel('f (kHz)'), ylabel('dB'), 
axis([0, 6, -90, 10]),
grid on, title('power spectrum - noisy heartbeat')

%
% 2.3.b Design an FIR LPF to filter out noise We see that the dominant
% frequency of the noisy and filtered signals is at about 500 Hz, and
% there's a dip just after that, so we'll try a passband-edge frequency of
% 500 Hz. This worked with M = 1000, but after reducing the filter order
% the pass-band edge frequency needed to be lowered. Acceptable results
% were obtained with it set to 250 Hz and M = 200. Tried it with M = 100
% and M = 125, but the results were too noisy. M = 150 worked fine enough.
% Results closer to that in the figure given can be obtained by increasing
% the pass-band edge and increasing the filter order, but since the design
% specs are to minimize M, I chose to keep it at Fp = 250 Hz and M = 150.
% So the cut-off is not as steep, but the results are still acceptable and
% the order is lower.

N = length(x);
M   = 150;        % FIR filter order
Fp  = 250;        % 600 Hz passband-edge frequency
Fs  = fs;         % 22.05 kHz sampling frequency
Rp  = 0.00057565; % Corresponds to 0.01 dB peak-to-peak ripple
Rst = 1e-4;       % Corresponds to 80 dB stopband attenuation

Num = firceqrip(M,Fp/(Fs/2),[Rp Rst],'passedge'); % eqnum = vec of coeffs
% fvtool(Num,'Fs',Fs,'Color','White') % Visualize filter

yFilt = filter(Num,1,x);


% Spectrum of y
[Py,f] = pwelch(yFilt,[],[],[],fs);     % power spectrum estimate
Pdb = 10*log10(abs(Py));            % convert to dB

figure          % Note: There seems to be an offset of about -45 dB
plot(f/1000,Pdb,'r'); xlabel('f (kHz)'), ylabel('dB'), 
axis([0, 6, -135, -35]),
grid on, title('power spectrum - filtered heartbeat')

figure
plot(t,yFilt, 'r'), axis ([0, 6.5, -0.5, 1]), title('filtered heartbeat'),
xlabel('t (sec)')

yFilt = yFilt/max(abs(yFilt)); % normalize for audio output

audiowrite('Filtered_Heartbeat.wav',yFilt,fs);

%
% soundsc(yFilt,fs)

%
[x,fs] = audioread('M1+B.wav'); % Mother + baby
[y,fs] = audioread('M2.wav');   % Mother's thorax

N = length(x);
% fs = 22.05 [kHz]

[x, xMean] = zmean(x'); % Make x zero-mean and the mean value
[y, yMean] = zmean(y'); % Make y zero-mean and the mean value 

t = (1:length(x))/fs;       % time index scaled for plotting

figure
plot(t,x),title('mother + baby''s heartbeat sound'),
axis([0, 6.5, -0.4, 0.4]), xlabel('t (sec)'), grid on 

figure
plot(t,y),title('mother''s heartbeat sound'),
axis([0, 6.5, -0.4, 0.4]), xlabel('t (sec)'), grid on 

% Run LMS algorithm on x and y
N = length(x);
M = 500;              
mu = 1e-3;
h = zeros(M+1,1);   % initialize empty filter values
w = zeros(M+1,1);   % initialize internal state buffer
% LMS algorithm
for n = 1:N
    w(1) = y(n);    % load internal state buffer
    xHat(n) = h'*w; % use h and y(n) to estimate of x(n)
    e(n) = x(n) - xHat(n);  % error output
    h = h + 2*mu*e(n)*w;    % update filter weights
    w = [w(1); w(1:end-1)]; % update buffer
end

figure
plot(t,e,'r'), title('baby''s heartbeat sound'),xlabel('t (sec)')
axis([0, 6.5, -0.4, 0.4]), xlabel('t (sec)'), grid on 

e = e/max(abs(e));

audiowrite('LMS_Heartbeat.wav',e,fs);

%
% soundsc(e,fs)