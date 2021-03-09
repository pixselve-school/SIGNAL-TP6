function [DSPM, f] = DensSpecPuiss(br, fs)
%Calcul la densité spectrale de puissance
%br est le bruit et fs la frequence d'echantillonnage
DSPM = size(br,1)*pwelch(br,1024); 
DSPM=DSPM(1:size(DSPM)/2);
L = size(DSPM,1);
f = fs*(1:L)/L-fs/2;
end

