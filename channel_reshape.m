function H_MIMO = channel_reshape(H, P)

Ntx = P.Ntx;
Nrx = P.Nrx;
CL = P.ChannelLength;
n_users = P.CDMAUsers;

H_MIMO = zeros(Nrx*CL, Ntx, n_users);
for usr = 1:n_users
    for c = 1:CL
        H_MIMO(1+(c-1)*Nrx:c*Nrx, :, usr) = H(:, : ,c , usr);
    end
end

end