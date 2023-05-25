function s_hat = MMSE_Detector_Biased(H, y, Pn)
% MMSE Detector
%   - for flat channel H and receive vector y
%   - slicing / quantization to constellations
%   - computes MMSE Moore-Penrose pseudo-inverse


N_tx = size(H,2);                       % # TX
N_rx = size(H,1);                       % # RX --> computed to account for cases in which N_tx ≠ N_rx

% calculate MMSE Moore-Penrose pseudo-inverse
Hi = (H'*H + Pn*N_tx*eye(N_tx))\H';
%Hi = diag(1./diag(Hi*H))*Hi;    % equalize power: diag(Hi*H) = [1 1...

% spatial MMSE equalization
y_hat = (Hi*y).';

% slicing / quantization to constellations
s_hat = real(y_hat(:));

end
