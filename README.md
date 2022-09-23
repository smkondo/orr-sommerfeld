# orrsommerfeld - laminar flow
Solving the Orr-sommerfeld problem as a generalized eigenvalue problem 

In this repository, the code to compute the eigenspectra and eigenfunction of the Orr-Sommerfeld problem are given. Please download all codes and run the "OSeigproblem.m" file to compute the eigenspectra and eigenfunction of the Orr-Sommerfeld problem.

# turbulent channel flow
The following MATLAB files represent the functions to run a resolvent framework for a Turbulent Channel Flow. The resolvent framework takes a mean streamwise velocity profile of a channel flow and returns the most energetic velocity fluctuations within the flow. This is an operator based modal analysis technique for turbulent flow calculations. 

A simple test case for Re_tau = 180 can be run by downloading all .m files and running test.m. In this test case, a first-rank approximatino of each fluctuating velocity components are computed and displayed (The figures can be found in the results section of this markdown file). This test case can easily be manipulated to use a different mean velocity profile, change the number of singular values, etc. 

## General Description of each .m file
A general description of each .m file is as follows:
* cheb_basis: value of Chebyshev polynomial of order n evaluated at each Chebyshev point y
* cheb_expansion_soln: Chebyshev expansion given the Chebyshev coefficients from SVD of the resolvent operator
* Der: Chebyshev collocation derivative matrix
* Deven: Computes matrix that converts Chebyshev coefficients of a polynomial to coefficients of derivative
* Dmat: Derivative matrices of the Chebyshev polynomials
* Fourier2physical: Converts from Fourier space to physical space
* meanU: computes an empircal streamwise mean velocity profile using the van Driest damping function
* meanUDNS: computes the streamwise mean velocity profile at each Chebyshev point using DNS mean velocity profile
* pois2: Computes the mass matrix (M) and linear operator (L) for computing the Orr-Sommerfeld and Squire operators
* resolvent: computes the first-rank approximation of the Fourier coefficients of the velocity fluctuations
* test: Plots the first-rank approximation of the fluctuating veloity contours at Re_tau = 1000 
* two: L2-norm weight of the Chebyshev polynomials

## Descrition of each figure 
1. Comparison between empirical and DNS mean velocity profile (Implementing empirical mean velocity profile found in Eqn. 3.2 from "Model-based scaling of the streamwise energy density in high-Reynolds number turbulent channels" by Rashad Moarref, Ati S. Sharma, Joel A. Tropp, and Beverley J. McKeon (2013))
![alt text](https://github.com/smkondo/orr-sommerfeld/blob/main/figures/fig1.png "mean velocity profile")
... The empirical formula for the mean velocity profile was derived by Reynolds and Tiederman (1967). The DNS data was provided by Lee and Moser (2015) and can be found [here](https://turbulence.oden.utexas.edu/channel2015/data/LM_Channel_0180_mean_prof.dat)

2. The 20 largest singular values (Recreating Fig. 4a from Moarref et. al. (2013))
![alt text](https://github.com/smkondo/orr-sommerfeld/blob/main/figures/fig2.png "singular values")
... As shown in the figure, the first rank approximation of the resolvent framework is valid since the first pair of singular values are approx. 6.3 times larger than the second pair of singular pairs. 

3. First-rank resolvent modes of fluctuating velocities (Recreating Fig. 7a from "Resolvent modelling of near-wall coherent structures in turbulent
channel flow" by Leandra I. Abreu, Andre V.G. Cavalieri, Philipp Schlatter, Ricardo Vinuesa, and Dan S. Henningson (2020))
![alt text](https://github.com/smkondo/orr-sommerfeld/blob/main/figures/fig3.png "resolvent modes")

4. First-rank approximation of streamwise velocity fluctuation (Recreating Fig. 5b from Abreu et. al. (2020))
![alt text](https://github.com/smkondo/orr-sommerfeld/blob/main/figures/fig4.png "streamwise fluctuation")

## References
The pdf files of the referenced literature can be found in this [google drive](https://drive.google.com/drive/folders/1EPNruFN7BVWhb943XOv2y26Jg5B-Fwzp?usp=sharing)