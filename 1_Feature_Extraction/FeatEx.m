% Delta	1
% Delta	2
% Theta	3
% Theta	4
% Alpha	5
% Alpha	6
% Sigma	7
% Sigma	8
% Beta	9
% Beta	10
% Gamma	11
% Delta	12
% Theta	13
% Alpha	14
% Sigma	15
% Beta	16
% Delta	17
% Theta	18
% Alpha	19
% Sigma	20
% Beta	21
% Delta	22
% Theta	23
% Alpha	24
% Sigma	25
% Beta	26
% Delta	35
% Delta	36
% Theta	37
% Theta	38
% Alpha	39
% Alpha	40
% Sigma	41
% Sigma	42
% Beta	43
% Beta	44
% Delta	45
% Delta	46
% Theta	47
% Theta	48
% Alpha	49
% Alpha	50
% Sigma	51
% Sigma	52
% Beta	53
% Beta	54




%% Feature Extraction


Fnum=69;       
Feats = zeros(1,Fnum);
n = length(Sig);  %number of samples in each epoch sample

X = Sig;
%% Power related features: RSP SWI and HP  total number 28
%================= RSP computation
Ptot=bandpower(X,Fs,[0,45]);
Feats(1) = bandpower(X,Fs,[.5,2])/Ptot;     %RSP_Delta1
Feats(2) = bandpower(X,Fs,[2,4])/Ptot;      %RSP_Delta2
Feats(3) = bandpower(X,Fs,[4,6])/Ptot;      %RSP_Theta1
Feats(4) = bandpower(X,Fs,[6,8])/Ptot;      %RSP_Theta2
Feats(5) = bandpower(X,Fs,[8,10])/Ptot;     %RSP_Alpha1
Feats(6) = bandpower(X,Fs,[10,12])/Ptot;    %RSP_Alpha2
Feats(7) = bandpower(X,Fs,[12,14])/Ptot;    %RSP_Sigma1
Feats(8) = bandpower(X,Fs,[14,16])/Ptot;    %RSP_Sigma2
Feats(9) = bandpower(X,Fs,[16,24])/Ptot;    %RSP_Beta1
Feats(10)= bandpower(X,Fs,[24,32])/Ptot;    %RSP_Beta2
Feats(11)= bandpower(X,Fs,[32,45])/Ptot;    %RSP_Gamma

%================= HP computation  Fc and S(Fc)
[P,f] = pwelch(X,[],[],0:10^-1:Fs/2,Fs);

%----> Fc    computation
Inds = (f>=0.5 & f<4);%Delta
Feats(12)=sum(P(Inds).*f(Inds))/sum(P(Inds));               %Delta
Inds = (f>=4 & f<8);%Theta
Feats(13)=sum(P(Inds).*f(Inds))/sum(P(Inds));               %Theta
Inds = (f>=8 & f<12);%Alpha
Feats(14)=sum(P(Inds).*f(Inds))/sum(P(Inds));               %Alpha
Inds = (f>=12 & f<16);%Sigma
Feats(15)=sum(P(Inds).*f(Inds))/sum(P(Inds));               %Sigma
Inds = (f>=16 & f<32);%Beta
Feats(16)=sum(P(Inds).*f(Inds))/sum(P(Inds));               %Beta
Inds = (f>=32 & f<45);%Gamma
Feats(17)=sum(P(Inds).*f(Inds))/sum(P(Inds));               %Gamma

%----> Fsigma    computation
Inds = (f>=0.5 & f<4);%Delta
Feats(18)=sqrt(sum(P(Inds).*(f(Inds)-Feats(12)).^2)/sum(P(Inds)));
Inds = (f>=4 & f<8);%Theta
Feats(19)=sqrt(sum(P(Inds).*(f(Inds)-Feats(13)).^2)/sum(P(Inds)));
Inds = (f>=8 & f<12);%Alpha
Feats(20)=sqrt(sum(P(Inds).*(f(Inds)-Feats(14)).^2)/sum(P(Inds)));
Inds = (f>=12 & f<16);%Sigma
Feats(21)=sqrt(sum(P(Inds).*(f(Inds)-Feats(15)).^2)/sum(P(Inds)));
Inds = (f>=16 & f<32);%Beta
Feats(22)=sqrt(sum(P(Inds).*(f(Inds)-Feats(16)).^2)/sum(P(Inds)));
Inds = (f>=32 & f<45);%Gamma
Feats(23)=sqrt(sum(P(Inds).*(f(Inds)-Feats(17)).^2)/sum(P(Inds)));

%----> S(fc) computation
Feats(24)=P(round(10*Feats(12))+1);
Feats(25)=P(round(10*Feats(13))+1);
Feats(26)=P(round(10*Feats(14))+1);
Feats(27)=P(round(10*Feats(15))+1);
Feats(28)=P(round(10*Feats(16))+1);
Feats(29)=P(round(10*Feats(17))+1);

%================= SWI computation
bspD = bandpower(X,Fs,[.6,4]);% subband spectrul power of delta
bspD2 = bandpower(X,Fs,[2,4]);% subband spectrul power of delta
bspT = bandpower(X,Fs,[4,8]);% subband spectrul power of theta
bspA = bandpower(X,Fs,[8,11.5]);% subband spectrul power of aplha based on Jobert et. al. 1994
%----> DSI    computation
Feats(30)=bspD/(bspT+bspA);
%----> TSI    computation
Feats(31)=bspT/(bspD+bspA);
%----> ASI    computation
Feats(32)=bspA/(bspD2+bspT);
%% Hjorth features
Xp=diff(X)*Fs;         %1st Differentiate
Xpp=diff(X,2)*Fs^2;    %2nd Differentiate
%----> Activity
Feats(33)=var(X);
%----> Mobility
Feats(34)=sqrt(var(Xp)/var(X));
%----> Complexity
Feats(35)=sqrt(var(Xpp)*var(X))/var(Xp);


%% Skewness and Kortusis
M2=0; %Moment 2
M3=0; %Moment 3  
M4=0; %Moment 4
m=mean(X);
for i=1:n
    M2=M2+(X(i)-m)^2;
    M3=M3+(X(i)-m)^3;
    M4=M4+(X(i)-m)^4;
end
M2=M2/n;
M3=M3/n;
M4=M4/n;
%----> Skewness
Feats(36)=M3/sqrt(M2^3);
%----> Kurtosis
Feats(37)=M4/M2^2;     

%% Bispectrum Time Series
[BTS,FFF] = btsestimate(X',Phi,Rho,nLag,Fs,nFFT);
desiredFreqs = [1,3,5,7,9,11,13,15,20,28,36,40];
[~,IndsF] = findpeaks(-min(abs(repmat(FFF,length(desiredFreqs),1)-repmat(desiredFreqs',1,length(FFF)))));
Feats(38:49) = abs(BTS(IndsF));
Feats(50:61) = angle(BTS(IndsF));

%% Wavelet Coefficients
WaveletLvl = 8;
Wavename = 'sym6';
type = 'd';
[C,L] = wavedec(X,WaveletLvl,Wavename);
PwaveletCoefs = zeros(1,WaveletLvl);
for i=1:WaveletLvl
    PwaveletCoefs(i) = bandpower(wrcoef('d',C,L,Wavename,i))/Ptot;
end
Feats(62:69)=PwaveletCoefs;

