close all;clear;
%% parameter initial
DSMorder = 4;
DSM_Type = 1;

BitPerSym = 10;   % OFDM format
FiberLen = 0e3*1;    % fiber length
Up = 8;         % OSR of delta-sigma modulator  
Up2 = 1;       % for PAM
SignalBW = 4e9;
SampleRateDefault = SignalBW*Up;
quant = 0;

OSNR_all = [50];
SNR_all = OSNR_all +10*log10(12.5e9/SampleRateDefault);
BER_all = [];
constellation = qammod([0:2^BitPerSym-1].',2^BitPerSym,'InputType','Integer').';
P = sqrt(mean(abs(constellation).^2));
%%  ---------- transmitter------------------ %%
%% data generation & OFDM modulation
N_fft = 1000*Up;    % fft point in each OFDM symbol
N_eff = 450;          % effective subcarrier number
N_zero = 7;           % subcarrier number around zero frequency
N_grd = N_fft - N_eff*2 - N_zero; 
OFDM_Num = 80;  % number of OFDM symbol 
rand('seed', 12);
TxBit = randi([0 1],1, N_eff*OFDM_Num*BitPerSym);
TxSym = qammod(TxBit.',2^BitPerSym,'InputType','Bit').';
TxSym_Mat0 = reshape(TxSym, N_eff, OFDM_Num);
TxSym_Mat = [zeros(floor((N_zero+1)/2),OFDM_Num);TxSym_Mat0;zeros(N_grd,OFDM_Num);conj(flipud(TxSym_Mat0));zeros(floor((N_zero-1)/2),OFDM_Num)];
TxWfm = ifft(TxSym_Mat);
TxWfm = reshape(TxWfm, 1, []);

TxWfm = TxWfm./max(abs(TxWfm));

PAPR = 10*log10(max(abs(TxWfm).^2)/mean(abs(TxWfm).^2))
figure; plot(10*log10(abs(fftshift(fft(TxWfm)))));   % modulated signal spectrum
%% DSM
switch DSM_Type
    case 0
        switch DSMorder
            case 1
                [ TxWfm ] = DSM1( TxWfm, 0.4, 2, 2 );
            case 2
                [ TxWfm ] = DSM2( TxWfm, 0.4, 2, 2 );
            case 3
                [ TxWfm ] = DSM3_v0( TxWfm, 0.4, 2, 2);
            case 4
                [ TxWfm ] = DSM4_CRFF( TxWfm, 0.4, 2, 2);
        end
        
    case 1
        AmplitudeRatio0 = 1.15*1;
        AmplitudeRatio1 = 0.2*1;
        [ TxWfm00 ] = DSM4_CRFF( TxWfm, AmplitudeRatio0, 2, 1);
        TxWfmFFT = fft(TxWfm00 - TxWfm);
        LLen = length(TxWfmFFT);
        kkk = 1.1;
        TxWfmFFT(round(LLen/Up/2*kkk)+2:end-round(LLen/Up/2*kkk)) = 0;
        TxWfmDel = ifft(TxWfmFFT);
        [ TxWfm01 ] = DSM4_CRFF( TxWfmDel, AmplitudeRatio1, 2, 1);
        
        TxWfm = TxWfm00*2 + TxWfm01; 
end
%% add preamble
rand('seed', 1);
TrainLen = 1024;
TrainBit = randi([0 2^2-1],1, TrainLen);
TrainSeq = pammod(TrainBit.',2^2,0,'gray').';
SynLen = 128;
SynSeq = TrainSeq(1:SynLen);
TxWfm = [TrainSeq TxWfm];
%% Mach-Zehnder modulator
TxWfm = TxWfm/sqrt(mean(abs(TxWfm).^2))*0.1;  % modulation index
Vpi = 6;                  % Vpi
Vbias = 4.5;    % bias
TxWfm1 = exp(1i*(TxWfm*1+Vbias*1)/Vpi*pi) + exp(-1i*(TxWfm*1+Vbias*1)/Vpi*pi);  % MZM emulation

figure; plot(10*log10(abs(fftshift(fft(TxWfm1)))));   % modulated signal spectrum
%% Fiber dispersion
C_speed=3e8;          % speed of light
lamda = C_speed/193.1e12;    % wavelength
CD_value = 17e-6 * FiberLen;    % accumulated dispersion
N = length(TxWfm1);
TimeW = N/SampleRateDefault;
beta2 = -lamda^2*CD_value/2/pi/C_speed;
w = 2*pi*[(0:N/2-1),(-N/2:-1)]/TimeW;
TxWfm1_FFT = fft(TxWfm1);
TxWfm1_FFT = TxWfm1_FFT.*exp(1i*-beta2*(w.^2)/2);
TxWfm1 = ifft(TxWfm1_FFT);

%% AWGN
for jj = 1:length(SNR_all)
    close all;
    SNR = 10^(SNR_all(jj)/10);
    Psig = mean(abs(TxWfm1).^2);
    randn('seed',1);
    noise = ( randn(1,length(TxWfm1))+1i*randn(1,length(TxWfm1)) )*sqrt(Psig/SNR/2);
    RxWfm1 = TxWfm1 + noise*1;
%% ---------- receiver------------------ %%
%% PD detection
RxWfm2 = abs(RxWfm1).^2;
%% DC block
RxWfm2 = RxWfm2 - mean(RxWfm2);
RxSig1 = RxWfm2./sqrt(mean(abs(RxWfm2).^2));
figure; plot(linspace(-SampleRateDefault/2/1e9,SampleRateDefault/2/1e9,length(RxWfm2)),10*log10(abs(fftshift(fft(RxWfm2)))));
%% synchronization
corr1 = zeros(1, 40000);
for ii = 1:40000
    corr1(ii) = abs(RxSig1(ii : 1 : ii+1*(SynLen-1) ) * SynSeq');
end
figure; plot(corr1);
[value, pos] = max(corr1);
RxSyn1 = RxSig1(pos : 1 : end); 
%% Equalization
delay = 5;
lambda = 0.9999;
EqMode = 1;
switch EqMode
    case 0
        RxEqu = RxSyn1;
    case 1
        [ RxEqu ] = RLS_Train( RxSyn1, TrainSeq, delay, lambda, 1 );
end
RxEqu = RxEqu(TrainLen+1:end);
%%
RxEqu = PAMdecision(RxEqu,[-3 -1 1 3]);

switch DSM_Type
    case 0
%% Low-pass filter
RxWfm2 = RxEqu;
RxWfm2FFT = fft(RxWfm2);
cut = round(length(RxWfm2)/Up*1.1);
RxWfm2FFT(cut+2:end-cut) = 0;
RxWfm2 = ifft(RxWfm2FFT);
figure; plot(linspace(-SampleRateDefault/2/1e9,SampleRateDefault/2/1e9,length(RxWfm2)),10*log10(abs(fftshift(fft(RxWfm2)))));
    case 1
%% modified
RxWfm2 = RxEqu;
RxWfm20 = (RxWfm2>0)*2-1;
RxWfm21 = RxWfm2 - 2*RxWfm20;
        
RxWfm2 = RxWfm20*AmplitudeRatio0 - RxWfm21*AmplitudeRatio1*1;        

RxWfm2FFT = fft(RxWfm2);
cut = round(length(RxWfm2)/Up*0.65);
RxWfm2FFT(cut+2:end-cut) = 0;
RxWfm2 = ifft(RxWfm2FFT);
end
%% OFDM demodulation
RxWfm2 = RxWfm2(1:N_fft*OFDM_Num);
RxWfm_Mat = reshape(RxWfm2,N_fft,[]);
RxWfm_Mat = fft(RxWfm_Mat);
RxSym_Mat = RxWfm_Mat( floor((N_zero+1)/2)+1 : floor((N_zero+1)/2)+N_eff , : );
RxSym = reshape(RxSym_Mat,1,N_eff*OFDM_Num);
RxSym = RxSym./sqrt(mean(abs(RxSym).^2))*P;
%%
scatterplot(RxSym);
RxBit = qamdemod(RxSym.',2^BitPerSym,'OutputType','Bit').';
%% Cal BER
[NumErr, BER, individual] = biterr(RxBit(1:end),TxBit(1:end));
BER
BER_all = [BER_all BER];
%% Cal SNR
RxSym_Mat0 = reshape(RxSym,N_eff,OFDM_Num);
SNR_color = 10*log10( mean(abs(TxSym_Mat0).^2, 2) ./mean( abs(RxSym_Mat0 - TxSym_Mat0).^2, 2) );
% figure; plot(SNR_color);
% figure;plot(10.^(-SNR_color/20),'r');
SNR = 10*log10( mean(abs(TxSym).^2 ) ./mean( abs(RxSym - TxSym).^2) );
SNR
EVM = 10^(-SNR/20)
end
