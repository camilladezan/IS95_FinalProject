% Advanced Wireless Receivers - Final Project:
% 
% CDMA implementation
%
% Camilla De Zan, Linda Fabiani
%
% Telecommunications Circuits Laboratory
% EPFL
clear all, close all, clc

%% Part 1 - IS95-CDMA implementation
% Parameters
P.NumberOfFrames      = 1000;
P.NumberOfSymbols     = 172;

P.AccessType = 'CDMA'; 
P.CDMAUsers     = 2;

P.Modulation    = 1;                   % 1: BPSK

P.ChannelType   = 'AWGN';         % 'AWGN', 'Multipath', 'PassThrough'
P.ChannelLength = 1;                   % > 1 for multipath
P.CoherenceTime = 100;             

P.KConvDecoder = 9;     % parameter K for the Viterbi decoder, linked to the traceback depth
P.ConvEncRate = 2;      % inverse of the rate of the convolutional encoder
P.poly = [753 561];     % polynomial for Viterbi decoder
P.HamLen = 64;          % Length of Hadamard Sequence

P.SNRRange = -30:2:10; % SNR Range to simulate in dB

P.ReceiverType  = 'Rake';
P.RakeFingers = 2;

BER = simulator(P);

simlab = sprintf('%s - Length: %d - Users: %d' ,P.ChannelType,P.ChannelLength,P.CDMAUsers);


figure(1)
semilogy(P.SNRRange,BER,'b.-','DisplayName',simlab, 'LineWidth', 1.2)

xlabel('SNR','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
grid minor;
legend('-DynamicLegend');

%hold on 

%% Section to run to make comparisons plots

P.ChannelType = 'AWGN';
P.ChannelLength = 1;

BER2 = simulator(P);
simlab = sprintf('%s - Users: %d' ,P.ChannelType,P.CDMAUsers);
% add BER curve to figure
figure(2)
semilogy(P.SNRRange,BER2,'b.-','DisplayName',simlab, 'LineWidth', 1.2)

xlabel('SNR','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
grid minor;
legend('-DynamicLegend');

% P.RakeFingers = 9;
% BER3 = simulator(P);
% simlab = sprintf('%s - Length: %d - Fingers: %d' ,P.ChannelType,P.ChannelLength,P.RakeFingers);
% % add BER curve to figure
% semilogy(P.SNRRange,BER3,'g.-','DisplayName',simlab, 'LineWidth', 1.2)
