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
P.CDMAUsers     = 1;

P.Modulation    = 1;                   % 1: BPSK

P.ChannelType   = 'Multipath';         % 'AWGN', 'Fading', 'Multipath', 'PassThrough'
P.ChannelLength = 2;                   % > 1 for multipath
P.CoherenceTime = 100;             

P.KConvDecoder = 9;     % parameter K for the Viterbi decoder, linked to the traceback depth
P.ConvEncRate = 2;      % inverse of the rate of the convolutional encoder
P.poly = [753 561];     % polynomial for Viterbi decoder
P.HamLen = 64;          % Length of Hadamard Sequence

P.SNRRange = -30:2:10; % SNR Range to simulate in dB

P.ReceiverType  = 'Rake';
P.RakeFingers = 3;

BER = simulator(P);

simlab = sprintf('%s - Length: %d - Users: %d' ,P.ChannelType,P.ChannelLength,P.CDMAUsers);


figure(1)
semilogy(P.SNRRange,BER,'b.-','DisplayName',simlab)

xlabel('SNR','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
grid minor;
legend('-DynamicLegend');

hold on 

%% Section to run to make comparisons plots

P.ChannelLength = 6;
BER2 = simulator(P);
simlab = sprintf('%s - Length: %d - Users: %d' ,P.ChannelType,P.ChannelLength,P.CDMAUsers);
% add BER curve to figure
semilogy(P.SNRRange,BER2,'r.-','DisplayName',simlab)

P.ChannelLength = 9;
BER3 = simulator(P);
simlab = sprintf('%s - Length: %d - Users: %d' ,P.ChannelType,P.ChannelLength,P.CDMAUsers);
% add BER curve to figure
semilogy(P.SNRRange,BER3,'g.-','DisplayName',simlab)
