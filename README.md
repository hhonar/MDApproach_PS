# MDApproach_PS

This repository is focused on introducing the data-driven mode decomposition based phase synchronization framework as a way to investigate time-varying phase synchrony.  This allows for taking an in-depth look over this technique for assessing time-varying functional brain connectivity, or additionally any potential phase syncrhony between multivariate signals.


#### Highlights:
• Interest in measuring phase synchronization between brain regions across time.\
• Analysis necessitates signals be sufficiently narrow-bandpassed.\
• Using data-driven mode decomposition (MD) approaches to achieve this goal.\
• To compare various MD approaches.\
• Multivariate variational mode decomposition outperforms other MD approaches.

##### Description:
This repository, while is intended to be used as a supplement to the Ref [i] below for assessing the time varying functional brain connectivity, it can also be leveraged as a signal processing and data driven technique for other applications as well where the synchrony between signals/time series are of interest.  

1. Simulations 1-6 are named Sim123, Sim4, Sim5, and Sim6
2. bemd.m: bivariate empirical mode decomposition
3. namemd.m: noise assisted multivariate empirical mode decomposition
4. mvmd.m: multivariate variational mode decompsoition 


The app version of this repository in MATLAB which is based on the scripts here will be coming soon. 

Please cite the following if you use any part of the script or code in your analysis and the appropriate reference listed below.
  [i]  Honari, H. Lindquist, M. A. (2022) Mode decomposition-based time-varying phase synchronization for fMRI. NeuroImag Journal.\
  [ii] Honari, H., Choe, A. S., & Lindquist, M. A. (2021). Evaluating phase synchronization methods in fMRI: A comparison study and new approaches. NeuroImage, 228, 117704.\
  [ii] Honari, H., Lindquist, M. A. (2021). Measuring time-varying connectivity using tapered Windowed Phase Synchronization. 27th Annual Meeting of the Organization for Human Brain Mapping, Vol. 2, pp. 58
  
The codes provided here are written in MATLAB.  The codes might require Signal Processing toolbox and some dependencies.  The proper references and dependencies for the codes are cited in the reference of the paper.  

References:

  [1] N. Rehman, H. Aftab, Multivariate Variational Mode Decomposition, arXiv:1907.04509, 2019. \
  [2] K. Dragomiretskiy, D. Zosso, Variational Mode Decomposition, IEEE Transactions on Signal Processing, vol. 62, pp. 531-544, 2014. \
  [3] Y. Zhang, P. Xu,P.Y. Li,K.Y.Duan,Y.X.Wen,Q.Yang,T.Zhang, and D.Z.Yao, Noise-assisted multivariate empirical mode decomposition for multichannel EMG signals[J]. Biomedical Engineering Online, 2017, 16(1):107.\
  [4] Y. Zhang, S. Su, P. Xu, D.Z Yao, Performance Evaluation of Noise-Assisted Multivariate Empirical Mode Decomposition and its Application to Multichannel EMG Signals, 39th Annual International Conference of the IEEE Engineering in Medicine and Biology Society, 2017, Jeju, Korea\
  [5]  Rehman and D. P. Mandic, "Multivariate Empirical Mode Decomposition", Proceedings of the Royal Society A, 2010\
  [6]  G. Rilling, P. Flandrin and P. Goncalves, "On Empirical Mode Decomposition and its Algorithms", Proc of the IEEE-EURASIP
       Workshop on Nonlinear Signal and Image Processing, NSIP-03, Grado (I), June 2003\
  [7]  N. E. Huang et al., "A confidence limit for the Empirical Mode Decomposition and Hilbert spectral analysis",
       Proceedings of the Royal Society A, Vol. 459, pp. 2317-2345, 2003
  
