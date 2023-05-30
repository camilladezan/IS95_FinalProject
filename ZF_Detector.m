function s_hat = ZF_Detector(H, y)
% Zero Forcing Detector
%   - for flat channel H and receive vector y
%   - slicing / quantization to constellations
%   - computes ZF Moore-Penrose pseudo-inverse

% calculate ZF Moore-Penrose pseudo-inverse
Hi = inv(H'*H)*H';

% spatial ZF equalization
y_hat = (Hi*y).';

% slicing / quantization to constellations
s_hat = real(y_hat(:));

end

