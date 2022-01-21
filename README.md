# Fair-Matrix-Design-for-Linear-Video-Coders

Master's research project (2018-2019), with the goal of designing fair encoding and decoding matrices for linear video coders (where all image transformations are linear). We use linear coding methods to avoid a sharp decrease in the image reconstruction fidelity when the receiver channel quality (SNR) is below the threshold for which the coding scheme was designed. We work in a multicast scenario (one source and multiple receivers) and minimize the greatest reconstruction error among all receiving channels, under a power constraint. 

<p align="center">
<img src="SIMO.png" class="centerImage" alt="drawing" width="600"/>
 </p>

This is particularly interesting if you send a message to multiple outputs and want to make sure that all channels can decode it with a satisfactory quality. Below, you can see the algorithm's solution makes sure no channel has too bad a decoding.

<p align="center">
<img src="min-Max solution.jpg" alt="drawing" width="600"/>
</p>

References : 
M. R. A. Khandaker and Y. Rong, "Precoding Design for MIMO Relay Multicasting," in IEEE Transactions on Wireless Communications, vol. 12, no. 7, pp. 3544-3555, July 2013, doi: 10.1109/TWC.2013.060413.121817.

Shuo Zheng. Accounting for channel constraints in joint source-channel video coding schemes. Signal and Image Processing. Université Paris-Saclay, 2019. English. ⟨NNT : 2019SACLT005⟩. ⟨tel-02050971⟩.
