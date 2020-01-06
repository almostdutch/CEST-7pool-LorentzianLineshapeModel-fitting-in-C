# Program to generate CEST maps by fitting a 7-pool Lorentzian-lineshape model to densely-sampled brain CEST-MRI data.

The program is dependent on the following Libraries:

1. Levmar (provided with this program). Levenberg-Marquardt nonlinear least squares algorithm. Available from http://users.ics.forth.gr/~lourakis/levmar/index.html#download

2. Intel Math Kernel Library (Intel MKL), which is free for Linux. For alternatives, see http://users.ics.forth.gr/~lourakis/levmar/index.html#download

3. MATLAB's I/O Library for reading/writing .mat files in C. An alternative (a bit slower than the native MATLAB library) is to use https://github.com/tbeu/matio

Test MATLAB dataset (CEST data 1 slice, 93 frequency offsets) is provided in /matlab/test-data/testDataIn.mat together with the expected results obtained from running this program: 

testDataOut_c2m.mat is the final output file\
results.fig and results.jpg are the examples of fitted CEST maps for each pool modeled. 

To compile and run the program (also for further clarification) see /src/COMPILATIONS_INSTRUCTIONS.txt
 
