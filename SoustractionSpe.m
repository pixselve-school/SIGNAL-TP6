function [den,den1,den2,hamx1,hamx2] = SoustractionSpecorrection(x1,hamming,DSPM)
alpha = 1.5 * 0.1;
beta = 0.005 * 0.1;
N=272*128; % = 136*256
%N = floor(length(x1)/2);
N = floor(length(x1)/256);
hamx1=zeros(N,1);
hamx2=zeros(N,1);
den1 = zeros(size(x1));
den2 = zeros(size(x1));
icplx=sqrt(-1);
delta=0.5;
for b = 1:(floor(length(x1)/256)-1) %135-20
    for i=1:256
        hamx1(256*b+i) = hamming(i)*x1(256*b+i);
    end
    TFhamx1{b}= abs(fft(hamx1((256*b+1:256*b+256)))).^2;
    ModS1{b} = (TFhamx1{b}).^delta-alpha.*DSPM.^delta;
    test1(256*b+1:256*b+256) = (ModS1{b}<beta*DSPM);
    ModS1{b} = ModS1{b}.*(ModS1{b}>beta*DSPM)+beta*DSPM.*(ModS1{b}<beta*DSPM);
%     ModS1{b} = sqrt(ModS1{b});
    phix1{b} = angle(fft(x1(256*b+1:256*b+256)));
    S1chap{b} = ModS1{b}.*exp(icplx*phix1{b});
    den1(256*b+1:256*b+256) = ifft(S1chap{b});
end

for c = 1:(floor(length(x1)/256)-2)
    for i=1:256
        hamx2(256*c+128+i) = hamming(i).*x1(128+256*c+i);
    end
    TFhamx2{c}= abs(fft(hamx2((256*c+128+1:256*c+128+256)))).^2;
    ModS2{c} = TFhamx2{c}.^delta-alpha.*DSPM.^delta;
    ModS2{c} = ModS1{c}.*(ModS2{c}>beta*DSPM)+beta*DSPM.*(ModS2{c}>beta*DSPM);
  %  ModS2{c} = sqrt(ModS2{c});
    phix2{c} = angle(fft(x1(256*c+128+1:256*c+128+256)));
    S2chap{c} = ModS2{c}.*exp(icplx*phix2{c});
    den2(256*c+1:256*c+256) = ifft(S2chap{c});
    test2(256*c+1:256*c+256) = (ModS2{c}<beta*DSPM);
end

m1=max(abs(fft(x1)));
den = 0.5*(den1+den2);
m2=max(abs(fft(den)));
den=((m2)/(m1))*(den);
end

