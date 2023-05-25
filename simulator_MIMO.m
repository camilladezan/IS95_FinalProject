function BER = simulator_MIMO(P)

if P.CDMAUsers > P.HamLen
    disp('WARNING: More user then sequences');
    BER = -1;
    return;
end
RX = P.CDMAUsers;

% define the convolutional encoder and decoder
convEnc = comm.ConvolutionalEncoder('TerminationMethod','Terminated', 'TrellisStructure', poly2trellis(P.KConvDecoder, P.poly));

L_dec = P.KConvDecoder*5 - 1;       % traceback depth
convDec = comm.ViterbiDecoder('TerminationMethod','Terminated', 'TrellisStructure', poly2trellis(P.KConvDecoder, P.poly), 'TracebackDepth', L_dec);

% Generate the spreading sequences
HadamardMatrix = hadamard(P.HamLen)/sqrt(P.HamLen);
SpreadSequence = HadamardMatrix;

SeqLen         = P.HamLen;

NumberOfBits   = P.NumberOfSymbols*P.Modulation; % per user

% definition of the Barker
NumberOfChips  = (P.NumberOfSymbols*P.Modulation + P.KConvDecoder -1)*SeqLen;           % per Frame
PNSequence     = genPNsequence(NumberOfChips);                                          % -(2*step(GS)-1);

% Channel
switch P.ChannelType
    case 'MIMO'
        NumberOfChipsRX   = NumberOfChips+P.ChannelLength-1;
    otherwise
        NumberOfChipsRX = NumberOfChips;
end

Results = zeros(1,length(P.SNRRange));

for ii = 1:P.NumberOfFrames

    ii

    bits = randi([0 1], NumberOfBits, RX); % Random Data

    % Encode data with convolutional encoded
    bits_encoded = zeros(P.ConvEncRate*(NumberOfBits+P.KConvDecoder-1), RX);
    for kk=1:RX
        bits_encoded(:, kk) = step(convEnc, bits(:, kk));
    end

    % Modulation
    switch P.Modulation % Modulate Symbols
        case 1 % BPSK
            symbols = -(2*bits_encoded - 1);
        otherwise
            disp('Modulation not supported')
    end

    % distribute symbols on users and antennas
    len = length(bits_encoded(:, 1))/P.Ntx;
    SymUsers = zeros(RX, len, P.Ntx);
    for tx_antenna=1:P.Ntx
        SymUsers(:, :, tx_antenna) = symbols(1+(tx_antenna-1)*len:tx_antenna*len, :)';
    end
    
    waveform = zeros(NumberOfChips, P.Ntx);
    for tx_antenna=1:P.Ntx
        % multiply hadamard on each antenna
        txsymbols_tmp= SpreadSequence(:,1:RX) * SymUsers(:, :, tx_antenna);
        
        % apply Barker code
        waveform(:, tx_antenna) = txsymbols_tmp(:).*PNSequence;
    end

    
    % Channel
    switch P.ChannelType
        case 'MIMO'
            H = sqrt(1/2)*randn(P.Nrx, P.Ntx, P.ChannelLength, RX) + sqrt(-1/2)*randn(P.Nrx, P.Ntx, P.ChannelLength, RX);
            H_P = permute(H, [3 1 2 4]);
            H_MIMO = reshape(H_P, [], size(H_P,3), size(H_P,4));
        otherwise
            disp('Channel not supported')
    end

    %%%

    % Simulation
    snoise = ( randn(NumberOfChipsRX, P.Nrx, RX) + 1i* randn(NumberOfChipsRX, P.Nrx, RX) );

    % SNR Range
    for ss = 1:length(P.SNRRange)
        SNRdb  = P.SNRRange(ss);
        SNRlin = 10^(SNRdb/10);
        noise  = 1/sqrt(2*SeqLen*SNRlin) * snoise;

        y = zeros(NumberOfChipsRX, P.Nrx, RX);
        % Channel
        switch P.ChannelType
            case 'MIMO'
                for i = 1:RX
                    for rx_antenna = 1:P.Nrx
                        for tx_antenna = 1:P.Ntx
                            y(:, rx_antenna, i) = y(:, rx_antenna, i) + conv(waveform(:, tx_antenna), squeeze(H(rx_antenna, tx_antenna, :, i)));
                        end
                    end
                end
                y = y + noise;
            otherwise
                disp('Channel not supported')
        end

        % Receiver
        switch P.ReceiverType
            case 'Rake'
                FrameLength = NumberOfChips;
                rxbits = zeros(P.NumberOfSymbols, RX);

                for usr = 1:P.CDMAUsers
                    for rx_antenna = 1:P.Nrx
                        rx_symb_finger = zeros(P.Nrx * min(P.ChannelLength, P.RakeFingers), (P.NumberOfSymbols*P.ConvEncRate+P.ConvEncRate*(P.KConvDecoder-1))/P.Ntx);
                        for finger = 1:min(P.RakeFingers,P.ChannelLength)
                            y_bark = y(finger:finger+FrameLength-1, rx_antenna, usr).*PNSequence;
                            rx_symb = reshape(y_bark,SeqLen,[]);
                            rx_symb_finger((rx_antenna-1)*min(P.ChannelLength, P.RakeFingers)+finger, :) = SpreadSequence(:,usr).'*rx_symb;
                        end
                    end

                    rx_symbols = rx_symb_finger;

                    switch P.MIMOdetector
                        case 'ZF'
                            H_user = squeeze(H_MIMO(:,:,usr));
                            s_hat = ZF_Detector(H_user, rx_symbols);
                        case 'SIC'
                            H_user = squeeze(H_MIMO(:,:,usr));
                            s_hat = SIC_Detector(H_user, rx_symbols, P.Constellation);
                        case 'MMSE'
                            % compute noise power
                            Pn = 1/SeqLen*SNRlin;
                            H_user = squeeze(H_MIMO(:,:,usr));
                            s_hat = MMSE_Detector_Biased(H_user, rx_symbols, Pn);
                            
                        otherwise
                            disp('Receiver not supported')
                    end
                    % convolutional decoder
                    rxbits_decoded = step(convDec, s_hat(:));
                    rxbits(:, usr) = rxbits_decoded(1:end-(P.KConvDecoder-1));
                end
            otherwise
                disp('Source Encoding not supported')
        end

        % BER count
        errors = sum(sum(rxbits ~= bits));
        Results(ss) = Results(ss) + errors;
    end
end

BER = Results/(P.Ntx*NumberOfBits*P.CDMAUsers*P.NumberOfFrames);

end


% Function to generate the PN sequence 
function seq = genPNsequence(len)

pnseq = comm.PNSequence('Polynomial',[42 35 33 31 27 26 25 22 21 19 18 17 16 10 7 6 5 3 2 1 0], ...
    'Mask',[1 1 0 0 0 1 1 0 0 0 randi([0,1], 1,32)], 'InitialConditions', ...
    [zeros(1,41) 1], 'SamplesPerFrame', len);

seq = pnseq();

end

