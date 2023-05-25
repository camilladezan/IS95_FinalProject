function sHat = SIC_Detector(H, y, constellations)
% SIC_DETECTOR_BIASED - SIC detector. Performes a SIC detection.
%    [sHat, complexity] = SIC_Detector(H, y, Pn_dB, Constellations)
%    
%    Arguments:
%      H:              (matrix)          Channel matrix
%      y:              (vector)          Receive vector
%    Return values:
%      sHat:           (integer vector)  Estimate of transmit vector
%     
% Author(s): YOU
% Copyright (c) 2012 TCL.


% noise and channel dimensions extraction
nTx = size(H,2);
n_data = length(y(1,:));

sHat = zeros(nTx, n_data);

Htmp = H;
ytmp = y;
for kk = 1:nTx
  % linear estimator
  G = (Htmp'*Htmp)\Htmp';
  sTilde = G(kk,:)*ytmp.';
    
  % symbol detection
  sHat(kk, :) = real(sTilde);

  % interference cancellation
  ytmp = ytmp - Htmp(:,kk)*sHat(kk, :);
  Htmp = Htmp(:,kk:end);
end


sHat = real(sHat(:));

end

