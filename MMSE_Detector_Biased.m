function s_hat = MMSE_Detector_Biased(H, y, Pn_dB, Constellations)
% MMSE Detector
%   - for flat channel H and receive vector y
%   - slicing / quantization to constellations
%   - computes MMSE Moore-Penrose pseudo-inverse

% compute noise power (relative to power 1)
Pn = 10^(Pn_dB/10);      % (Pn_dB is basically negative, but may be positive)
N_tx = size(H,2);                       % # TX
N_rx = size(H,1);                       % # RX --> computed to account for cases in which N_tx ≠ N_rx

% calculate MMSE Moore-Penrose pseudo-inverse
Hi = (H'*H + Pn*N_tx*eye(N_tx))\H';
%Hi = diag(1./diag(Hi*H))*Hi;    % equalize power: diag(Hi*H) = [1 1...

% spatial MMSE equalization
y_hat = Hi*y;

x_hat = zeros(size(y));
s_hat = zeros(size(H,2),size(y,2));

% slicing / quantization to constellations
for t = 1:size(y,2)
  for ntx = 1:N_tx
    tmp=Constellations-y_hat(ntx,t);
    tmp=abs(tmp).^2;
    [~,idx]=min(tmp);
    x_hat(ntx,t)=Constellations(idx);
    s_hat(ntx,t)=idx-1;
  end
end

end
