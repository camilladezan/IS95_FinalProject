% Advanced Wireless Receivers - Final Project:
%
% Camilla De Zan, Linda Fabiani
%
% Telecommunications Circuits Laboratory
% EPFL
clear all, close all, clc

% Parameters
P.NumberOfFrames      = 1000;
P.NumberOfSymbols     = 172;

P.AccessType = 'CDMA'; 
P.CDMAUsers     = 63;

P.Modulation    = 1;        % 1: BPSK

P.ChannelType   = 'Multipath';      % 'AWGN', 'Fading', 'Multipath', 'PassThrough'
P.ChannelLength = 3;                % increase it for multipath

P.KConvDecoder = 9;     % parameter K for the Viterbi decoder, linked to the traceback depth
P.ConvEncRate = 2;      % inverse of the rate of the convolutional encoder
P.poly = [753 561];     % polynomial for Viterbi decoder
P.HamLen = 64;          % Length of Hadamard Sequence

P.SNRRange = -15:2:10; % SNR Range to simulate in dB

P.ReceiverType  = 'Rake';
P.RakeFingers = 5;

BER = simulator(P); % simulator_coded_CDMA()/simulator_MIMO_coded_CDMA()....
% write simulators as matlab functions in separate files

simlab = sprintf('%s - Length: %d - Users: %d' ,P.ChannelType,P.ChannelLength,P.CDMAUsers);


hold on
semilogy(P.SNRRange,BER,'b.-','DisplayName',simlab)

xlabel('SNR','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
%xlim([min(P.SNRRange) max(P.SNRRange)]);
grid minor;
legend('-DynamicLegend');