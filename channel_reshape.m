function H_MIMO = channel_reshape(H, P)

Ntx = P.Ntx;
Nrx = P.Nrx;
CL = P.ChannelLength;
n_users = P.CDMAUsers;

H_MIMO = zeros(Nrx*CL, Ntx, n_users);
% for usr = 1:n_users
%     for c = 1:CL
%         H_MIMO(1+(c-1)*Nrx:c*Nrx, :, usr) = H(:, : ,c , usr);
%     end
% end

for usr = 1:n_users
    for tx = 1:P.Ntx
        for rx = 1:P.Nrx
            H_MIMO(1+(rx-1)*CL:rx*CL, tx, usr) = H(rx, tx, :, usr);
        end
    end
end

end