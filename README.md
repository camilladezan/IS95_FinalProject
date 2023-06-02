# IS95 FinalProject
Advanced Wireless Receivers Final Project: Implementation of the IS 95 standard for MIMO systems.

## Project description
Following the IS-95 standard we implement two simulators for the forward channel (i.e. from basestation to mobile) in MATLAB:
- the [first](https://github.com/camilladezan/IS95_FinalProject/blob/main/FinalProject.m) one deals with the single antenna scenario
- the [second](https://github.com/camilladezan/IS95_FinalProject/blob/main/FinalProject_MIMO.m) one extends the previous implementation to the MIMO case

The main building blocks of the implementation are a convolutional encoder and the corresponding Viterbi decoder, the CDMA spreading and despreading anf the RAKE receiver. In the MIMO implementation three types of detectors are tested: Zero-Forcing (ZF), Minimum Mean Squared Error (MMSE) and Successive Interference Cancellation (SIC-ZF).
The used parameters come from the standard documentation that can be found online. 

### Authors
Camilla De Zan and Linda Fabiani
