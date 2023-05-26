function sHat = SIC_Detector(H, y)
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
n_sym = length(y(1,:));

sHat = zeros(nTx, n_sym);

Htmp = H;
ytmp = y;
for kk = nTx:-1:1
  %G = (Htmp'*Htmp)\Htmp';
  G = pinv(Htmp);
  sTilde = G(kk, :)*ytmp;
  
  plot(sTilde, '.')
  % symbol detection
  sHat(kk, :) = real(sTilde);
  
  % interference cancellation
  ytmp = ytmp - Htmp(:, kk)*sHat(kk,:);
  Htmp = Htmp(:,1:kk-1);
end


sHat = real(sHat(:));

end

