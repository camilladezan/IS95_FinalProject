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
    
    Results = zeros(1,length(P.SNRRange));

    for ii = 1:P.NumberOfFrames
    
        %ii
    
        bits = randi([0 1],RX,NumberOfBits); % Random Data
    
        % Encode data with convolutional encoded
        for kk=1:RX
            bits_encoded(kk,:) = step(convEnc, bits(kk,:)');
        end
        
        % Modulation
        switch P.Modulation % Modulate Symbols
            case 1 % BPSK
                symbols = -(2*bits_encoded - 1);
            otherwise
                disp('Modulation not supported')
        end
       
        % distribute symbols on users and antennas
        % SymUsers = reshape(symbols, RX, [], P.Ntx);
        len = length(bits_encoded(1,:))/P.Ntx;
        SymUsers = zeros(RX, len, P.Ntx);
        for tx_antenna=1:P.Ntx
            SymUsers(:, :, tx_antenna) = symbols(:, 1+(tx_antenna-1)*len:tx_antenna*len);
        end
            
        % multiply hadamard on each antenna
        for tx_antenna=1:P.Ntx
            txsymbols(:, :, tx_antenna) = SpreadSequence(:,1:RX) * SymUsers(:, :, tx_antenna);
        end
    
        % definition of the Barker
        NumberOfChips  = length(txsymbols(:))/P.Ntx;          % per Frame
        PNSequence     = genbarker(NumberOfChips);            % -(2*step(GS)-1);

        % Channel
        switch P.ChannelType
            case 'MIMO'
                NumberOfChipsRX   = NumberOfChips+P.ChannelLength-1;
            otherwise
                NumberOfChipsRX = NumberOfChips;
        end

        
        % apply Barker code
        for tx_antenna=1:P.Ntx
            txsymbols_tmp = txsymbols(:, :, tx_antenna);
            waveform(tx_antenna, :) = txsymbols_tmp(:).*PNSequence;
        end
        
     
    
        % reshape to add multi RX antenna suppport
        waveform  = reshape(waveform, P.Ntx, NumberOfChips);
        mwaveform = repmat(waveform,[1 1 RX]);
        
        % Channel
        switch P.ChannelType
            case 'MIMO'
                H = sqrt(1/2)*randn(P.Nrx, P.Ntx, P.ChannelLength)+sqrt(-1/2)*randn(P.Nrx, P.Ntx, P.ChannelLength);
            otherwise
                disp('Channel not supported')
        end
    
        %%%
        
        % Simulation
        snoise = ( randn(P.Nrx,NumberOfChipsRX,RX) + 1i* randn(P.Nrx,NumberOfChipsRX,RX) );

        % SNR Range
        for ss = 1:length(P.SNRRange)
            SNRdb  = P.SNRRange(ss);
            SNRlin = 10^(SNRdb/10);
            noise  = 1/sqrt(2*SeqLen*SNRlin) *snoise;

            % Channel
            switch P.ChannelType
                case 'MIMO'
                    y = zeros(P.Nrx, NumberOfChipsRX, RX);
                    for rx_antenna = 1:P.Nrx
                        for tx_antenna = 1:P.Ntx
                            for i = 1:RX
                                y(rx_antenna,:,i) = y(rx_antenna,:,i) + conv(mwaveform(tx_antenna,:,i), squeeze(H(rx_antenna, tx_antenna, :)).');
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
                    FrameLength = NumberOfChipsRX;
                    rxbits = zeros(P.CDMAUsers,P.NumberOfSymbols);
                    for rx_antenna = 1:P.Nrx
                        for usr = 1:P.CDMAUsers
                            rx_symb_finger = zeros(P.Nrx, P.ChannelLength, (P.NumberOfSymbols*P.ConvEncRate+P.ConvEncRate*(P.KConvDecoder-1))/P.Ntx);
                            for finger = 1:min(P.RakeFingers,P.ChannelLength)
                                y_bark = y(rx_antenna,finger:finger+FrameLength-1,usr).*PNSequence.';
                                rx_symb = reshape(y_bark,SeqLen,[]);
                                rx_symb_finger(P.Nrx, finger, :) = SpreadSequence(:,usr).'*rx_symb;
                            end
                        end
                    end
                    % reshape the rx_symb_finger and the channel matrix to have an easier detector
                    % structure
                    rx_symbols = zeros(P.Nrx*P.ChannelLength, (P.NumberOfSymbols*P.ConvEncRate+P.ConvEncRate*(P.KConvDecoder-1))/P.Ntx);
                    for i = 1:min(P.RakeFingers,P.ChannelLength)
                        rx_symbols(1+(i-1)*P.Nrx:i*P.Nrx, :) = rx_symb_finger(:, i, :);
                    end
                    H_reshaped = zeros(P.Nrx*P.ChannelLength, P.Ntx);
                    for i = 1:P.Ntx
                        for j = 1:P.Nrx
                            H_reshaped(1+(j-1)*P.ChannelLength:j*P.ChannelLength, i) = H(j, i, :);
                        end
                    end
                    % scale constellation power
                    const_power=(norm(P.Constellation)^2)/2^P.Modulation;
                    scaling=sqrt(1/const_power);
                    P.Constellation=P.Constellation*scaling*sqrt(1/P.Ntx);
                    switch P.MIMOdetector 
                        case 'ZF'
                            s_hat = ZF_Detector(H, rx_symbols, P.Constellation);
                            s_hat = mexde2bi(s_hat(:), P.Modulation);
                        case 'SIC'
                            s_hat = SIC_Detector(H, rx_symbols, P.Constellation);
                            s_hat = mexde2bi(s_hat(:), P.Modulation);
                        case 'MMSE'
                            s_hat = MMSE_Detector(H, rx_symbols, P.Constellation);
                            s_hat = mexde2bi(s_hat(:), P.Modulation);
                        otherwise
                            disp('Receiver not supported')
                    end
                    % convolutional decoder
                    rxbits_encoded = step(convDec, s_hat);
                    rxbits(usr, :) = rxbits_encoded(1:end-(P.KConvDecoder-1));
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


function seq = genbarker(len)
    BarkerSeq = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];

    factor = ceil(len/length(BarkerSeq));
    b = repmat(BarkerSeq,1,factor);
    b = BarkerSeq.'*ones(1,factor);
    seq = b(1:len).';
end

