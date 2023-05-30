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

NumberOfBits   = P.NumberOfSymbols*P.Modulation;        % per user

Results = zeros(1,length(P.SNRRange));

for ii = 1:P.NumberOfFrames

    ii

    bits = randi([0 1], NumberOfBits, RX);      % Random Data generation

    % Encode data with convolutional encoded
    NumberOfEncodedBits = P.ConvEncRate*(NumberOfBits+P.KConvDecoder-1);
    bits_encoded = zeros(NumberOfEncodedBits, RX);
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

    % Create a matrix of symbols per user and spread symbols on TX antennas 
    SymUsers = zeros(NumberOfEncodedBits/P.Ntx, P.Ntx, RX);
    for usr = 1:RX
        SymUsers(:,:,usr) = reshape(symbols(:, usr), [], P.Ntx);
    end

    % multiply hadamard for the bits on each antenna
    txsymbols = zeros(P.HamLen, NumberOfEncodedBits/P.Ntx, P.Ntx);
    for tx_antenna = 1:P.Ntx
        txsymbols(:, :, tx_antenna) = SpreadSequence(:,1:RX) * SymUsers(:, tx_antenna, :).';
    end

    % definition of the PN sequence and BPSK modulation
    NumberOfChips  = size(txsymbols, 1)*size(txsymbols, 2);          % per Frame
    PNSequence     = -(2*genPNsequence(NumberOfChips)-1);      

    % Channel
    switch P.ChannelType
        case 'Multipath'
            NumberOfChipsRX   = NumberOfChips+P.ChannelLength-1;
        otherwise
            NumberOfChipsRX = NumberOfChips;
    end


    % apply PN sequence to each antenna
    waveform = zeros(NumberOfChips, P.Ntx);
    for tx_antenna = 1:P.Ntx
        txsymbols_tmp = txsymbols(:, :, tx_antenna);
        waveform(:, tx_antenna) = txsymbols_tmp(:).*PNSequence;
    end

    % reshape to duplicate the transmitted symbols for all the users in the
    % system
    mwaveform = repmat(waveform,[1 1 RX]);
     
    % Channel
    switch P.ChannelType
        case 'PassThrough'
            H = ones(P.Nrx, P.Ntx, RX);
            H_MIMO = H;
        case 'AWGN'
            H = ones(P.Nrx, P.Ntx, RX);
            H_MIMO = H;
        case 'Multipath'
            H = sqrt(1/2)*randn(P.Nrx, P.Ntx, P.ChannelLength, RX) + sqrt(-1/2)*randn(P.Nrx, P.Ntx, P.ChannelLength, RX);
            H_MIMO = channel_reshape(H, P);
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
            case 'PassThrough'
                for i =1:RX
                    y(:, :, i) = (H(:,:,i)*mwaveform(:,:,i).').';
                end
            case 'AWGN'
                for i =1:RX
                    y(:, :, i) = (H(:,:,i)*mwaveform(:,:,i).').';
                end
                y = y + noise;
            case 'Multipath'
                for i = 1:RX
                    for rx_antenna = 1:P.Nrx
                        for tx_antenna = 1:P.Ntx
                            y(:, rx_antenna, i) = y(:, rx_antenna, i) + conv(mwaveform(:, tx_antenna, i), squeeze(H(rx_antenna, tx_antenna, :, i)));
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
                rxbits = zeros(NumberOfBits, RX);
                rx_symb_finger = zeros(P.Nrx*min(P.ChannelLength, P.RakeFingers), FrameLength/P.HamLen);

                for usr = 1:P.CDMAUsers
                    for rx_antenna = 1:P.Nrx
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
                            s_hat = SIC_Detector(H_user, rx_symbols);
                        case 'MMSE'
                            % compute noise power
                            Pn = 1/(SeqLen*SNRlin);
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

    % Parameters
    tapPositions = [23 18 0];  % Tap positions for the feedback

    % Create PNSequence object
    pnGen = comm.PNSequence('Polynomial', tapPositions, 'InitialConditions', ones(1, 23), 'SamplesPerFrame', len);

    % Generate the PN sequence
    seq = pnGen();

end

