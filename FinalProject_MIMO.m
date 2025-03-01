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
P.NumberOfSymbols     = 172;

P.AccessType = 'CDMA'; 
P.CDMAUsers     = 2;

P.Modulation    = 1;        % 1: BPSK       
P.Constellation = [1 -1];

P.ChannelType   = 'Multipath';           % 'AWGN', 'Multipath', 'PassThrough'
P.ChannelLength = 3;                % increase it for multipath

P.KConvDecoder = 9;         % parameter K for the Viterbi decoder, linked to the traceback depth
P.ConvEncRate = 2;          % inverse of the rate of the convolutional encoder
P.poly = [753 561];         % polynomial for Viterbi decoder
P.HamLen = 64;              % Length of Hadamard Sequence

P.SNRRange = -30:2:10;      % SNR Range to simulate in dB

P.ReceiverType  = 'Rake';
P.RakeFingers = 4;

% definition of number of transmitting and receiving antennas
P.Ntx = 2; 
P.Nrx = 2;

% definition of the MIMO receiver
P.MIMOdetector = 'ZF';      %ZF, SIC, MMSE

BER = simulator_MIMO(P);

simlab = sprintf('%s - N_txN_r: %dx%d' ,P.MIMOdetector, P.Ntx, P.Nrx);

figure(1)
semilogy(P.SNRRange,BER,'b.-','DisplayName',simlab)
title('Detectors comparison')
xlabel('SNR','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
%xlim([min(P.SNRRange) max(P.SNRRange)]);
grid minor;
legend('-DynamicLegend');

%% Code to make comparison plots
hold on 

P.MIMOdetector = 'MMSE';
BER2 = simulator_MIMO(P);
simlab = sprintf('%s - N_txN_r: %dx%d' ,P.MIMOdetector, P.Ntx, P.Nrx);
semilogy(P.SNRRange,BER2,'r.-','DisplayName',simlab)


P.MIMOdetector = 'SIC';
BER3 = simulator_MIMO(P);
simlab = sprintf('%s - N_txN_r: %dx%d' ,P.MIMOdetector, P.Ntx, P.Nrx);
semilogy(P.SNRRange,BER3,'g.-','DisplayName',simlab)