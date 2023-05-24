function sHat = SIC_Detector(H, y, constellations)
% SIC_DETECTOR_BIASED - SIC detector. Performes a SIC detection.
%    [sHat, complexity] = SIC_Detector(H, y, Pn_dB, Constellations)
%    
%    Arguments:
%      H:              (matrix)          Channel matrix
%      y:              (vector)          Receive vector
%      Pn_dB:          (real)            Noise power
%      Constellations: (vector)          Vector with all possible constellations
%    Return values:
%      sHat:           (integer vector)  Estimate of transmit vector
%      complexity:     (integer)         Complexitx messure. (Not used)
% 
% Author(s): YOU
% Copyright (c) 2012 TCL.


% noise and channel dimensions extraction
nRx = size(H,1);
nTx = size(H,2);

sHat = zeros(nTx,1);
sIdx = zeros(nTx,1);

Htmp = H;
ytmp = y;
for kk = nTx:-1:1
  G = (Htmp'*Htmp)\Htmp';
  sTilde = G(kk,:)*ytmp;
  [~, constIdx] = min(abs(constellations - sTilde));
  sIdx(kk) = constIdx - 1;
  sHat(kk) = constellations(constIdx);
  ytmp = ytmp - Htmp(:,kk)*sHat(kk);
  Htmp = Htmp(:,1:kk-1);
end


sHat = sIdx;

end

