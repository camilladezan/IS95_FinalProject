function BER = simulator(P)

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
   
    % distribute symbols on users
    SymUsers = symbols;
        
    % multiply hadamard
    txsymbols = SpreadSequence(:,1:RX) * SymUsers;

    % definition of the Barker
    NumberOfChips  = length(txsymbols(:));          % per Frame
    PNSequence     = genbarker(NumberOfChips);      % -(2*step(GS)-1);

    % Channel
    switch P.ChannelType
        case 'Multipath'
            NumberOfChipsRX   = NumberOfChips+P.ChannelLength-1;
        otherwise
            NumberOfChipsRX = NumberOfChips;
    end

        
    % apply Barker code
    waveform = txsymbols(:).*PNSequence;

    % reshape to add multi RX antenna suppport
    waveform  = reshape(waveform,1,NumberOfChips);
    mwaveform = repmat(waveform,[1 1 RX]);
    
    % Channel
    switch P.ChannelType
        case 'PassThrough'
            himp = ones(RX,1);
        case 'AWGN'
            himp = ones(RX,1);
        case 'Multipath'
            himp = sqrt(1/2)* ( randn(RX,P.ChannelLength) + 1i * randn(RX,P.ChannelLength) );
        otherwise
            disp('Channel not supported')
    end
    
    %%%
    % Simulation
    snoise = ( randn(1,NumberOfChipsRX,RX) + 1i* randn(1,NumberOfChipsRX,RX) );
    
    % SNR Range
    for ss = 1:length(P.SNRRange)
        SNRdb  = P.SNRRange(ss);
        SNRlin = 10^(SNRdb/10);
        noise  = 1/sqrt(2*SeqLen*SNRlin) *snoise;
        
        % Channel
        switch P.ChannelType
            case 'PassThrough'
                y = mwaveform;
            case 'AWGN'
                y = mwaveform + noise;
            case 'Multipath'     
                y = zeros(1,NumberOfChips+P.ChannelLength-1,RX);
                for i = 1:RX
                    y(1,:,i) = conv(mwaveform(1,:,i),himp(i,:)) + noise(1,:,i); 
                end
            otherwise
                disp('Channel not supported')
        end
 
        % Receiver
        switch P.ReceiverType
            case 'Rake'
                FrameLength = NumberOfChips;
                rxbits = zeros(P.CDMAUsers,P.NumberOfSymbols);
                for usr = 1:P.CDMAUsers
                    rx_symb_finger = zeros(P.ChannelLength,P.NumberOfSymbols*P.ConvEncRate+P.ConvEncRate*(P.KConvDecoder-1));
                    for finger = 1:min(P.RakeFingers,P.ChannelLength)
                        y_bark = y(1,finger:finger+FrameLength-1,usr).*PNSequence.';
                        rx_symb = reshape(y_bark,SeqLen,[]);
                        rx_symb_finger(finger,:) = SpreadSequence(:,usr).'*rx_symb;
                    end
                    mrc = 1/norm(himp(usr,:))^2*conj(himp(usr,1:finger))*rx_symb_finger;
                    mrc = real(mrc);
                    % convolutional decoder
                    rxbits_encoded = step(convDec, mrc');
                    rxbits(usr, :) = rxbits_encoded(1:end-(P.KConvDecoder-1));

                end
            otherwise
            disp('Source Encoding not supported')
        end
        
        % BER count
        errors = sum(sum(rxbits ~= bits));
        Results(ss) = Results(ss) + errors;
    end
end

BER = Results/(NumberOfBits*P.CDMAUsers*P.NumberOfFrames);
end

function seq = genbarker(len)
    BarkerSeq = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];

    factor = ceil(len/length(BarkerSeq));
    b = repmat(BarkerSeq,1,factor);
    b = BarkerSeq.'*ones(1,factor);
    seq = b(1:len).';
end
