function s_hat = ZF_Detector(H, y, Constellations)
% Zero Forcing Detector
%   - for flat channel H and receive vector y
%   - slicing / quantization to constellations
%   - computes ZF Moore-Penrose pseudo-inverse

% calculate ZF Moore-Penrose pseudo-inverse
Hi = (H'*H)\H';

% spatial ZF equalization
y_hat = Hi*y;

x_hat = zeros(size(y));
s_hat = zeros(size(H,2),size(y,2));

% slicing / quantization to constellations
for t = 1:size(y,2)            % # data vector
  for ntx = 1:size(y_hat,1)    % # TX
    tmp=Constellations-y_hat(ntx,t);
    tmp=abs(tmp).^2;
    [~,idx]=min(tmp);
    x_hat(ntx,t)=Constellations(idx);
    s_hat(ntx,t)=idx-1;
  end
end

end

