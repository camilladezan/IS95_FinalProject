% Advanced Wireless Receivers - Final Project:
% 
% MIMO implementation
%
% Camilla De Zan, Linda Fabiani
%
% Telecommunications Circuits Laboratory
% EPFL

%% Part 2 - MIMO extension of the standard
close all, clear, clc

% Parameters
P.NumberOfFrames      = 1000;
P.NumberOfSymbols     = 125;

P.AccessType = 'CDMA'; 
P.CDMAUsers     = 1;

P.Modulation    = 1;        % 1: BPSK
P.Constellation = [1 -1];

P.ChannelType   = 'MIMO';           % 'AWGN', 'Fading', 'Multipath', 'PassThrough'
P.ChannelLength = 2;                % increase it for multipath

P.KConvDecoder = 9;         % parameter K for the Viterbi decoder, linked to the traceback depth
P.ConvEncRate = 2;          % inverse of the rate of the convolutional encoder
P.poly = [753 561];         % polynomial for Viterbi decoder
P.HamLen = 64;              % Length of Hadamard Sequence

P.SNRRange = -5:2:10;      % SNR Range to simulate in dB

P.ReceiverType  = 'Rake';
P.RakeFingers = 3;

% definition of number of transmitting and receiving antennas
P.Ntx = 2; 
P.Nrx = 2;

% definition of the MIMO receiver
P.MIMOdetector = 'ZF';      %ZF, SIC, MMSE

BER = simulator_MIMO(P);

simlab = sprintf('%s - Length: %d - Users: %d' ,P.ChannelType,P.ChannelLength,P.CDMAUsers);

figure(1)
semilogy(P.SNRRange,BER,'b.-','DisplayName',simlab)

xlabel('SNR','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
%xlim([min(P.SNRRange) max(P.SNRRange)]);
grid minor;
legend('-DynamicLegend');