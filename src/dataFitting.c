/*
 * Program to generate CEST maps by fitting a 7-pool Lorentzian-lineshape model
 *  	to densely-sampled brain CEST-MRI data.
 *
 * 7-pool Lorentzian-lineshape model:
 * (1) constant pool (offset). Accounts for difference in relaxation
 * 		properties between different tissues.
 * (2) MT pool. Magnetization transfer at -2.4 ppm.
 * (3) NOE1 pool. Nuclear Overhauser Enhancement at -3.5 ppm.
 * (4) NOE2 pool. Nuclear Overhauser Enhancement at -1.6 ppm.
 * (5) Cr pool. Cretine + other (amine) CEST effects at 2 ppm.
 * (6) Amide pool. (Poly)-peptides at 3.5 ppm.
 * (7) Water pool. Direct water saturation centered at 0 ppm.
 *
 * INPUT:
 *  The following options need to be specified for data fitting:
 * -i inputFileName.mat (MATLAB .mat file).
 * -d inputCestData (varible in inputFileName.mat).
 * 	inputCestData should have been normalized and B0-corrected.
 * -k originalFreqOffsets (varible in inputFileName.mat)
 * -o outputFileName.mat (MATLAB .mat file)
 * -l newFreqOffsets (varible in inputFileName.mat)
 * -s poolSwitchOnOff (varible in inputFileName.mat). A vector of 0 (exclude) and 1 (include)
 * 	for each of 7 pools in the model.
 * -r repetitions for data fitting
 * any (optional) comments (will be saved in file "outputFileName_logfile.txt")
 *
 * Example:
 * ./dataFitting -i inputFileName.mat -d inputCestData -k originalFreqOffsets
 *  	-o outputFileName.mat -l newFreqOffsets -s poolSwitchOnOff -r 100 my optional comments
 *
 * OUTPUT:
 * outputFileName.mat
 * Use c2matlab_data_conversion.m to convert data
 *  	from C style (row major order) to MATLAB style (column major order).
 * outputFileName_c2m.mat contains:
 * fitPars - Contains all fitted parameters [dim1 x dim2 x dim3 x 19], where
 * 1st parameter - constant offset
 * 2-4 parameters -   MT (amplitude, FWHM, offset)
 * 5-7 parameters -   NOE1 (amplitude, FWHM, offset)
 * 8-10 parameters -  NOE2 (amplitude, FWHM, offset)
 * 11-13 parameters - Cr (amplitude, FWHM, offset)
 * 14-16 parameters - Amide (amplitude, FWHM, offset)
 * 17-19 parameters - Water (amplitude, FWHM, offset)
 *
 * Individual fitted spectra for each pool are saved in the corresponding arrays:
 * fitMtPool -    MT
 * fitNoe1Pool -  NOE1
 * fitNoe2Pool -  NOE2
 * fitCrPool -    Cr
 * fitAmidePool - Amide
 * fitWaterPool - Water
 *
 * Time stamp and (optional) comments will be saved in file "outputFileName_logfile.txt"
 */

/*
 * (c) 2019, Vitaliy Khlebnikov, PhD
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <signal.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <unistd.h>
#include "levmar.h"
#include "misc.h"
#include "mat.h"
#include "matrix.h"
#include "compiler.h"
#include "multiLorentzian.h"

int main(int argc, char *argv[]) {

	if (catchSignal(SIGINT, interrupt) == -1) {
		error("Can't map the handler");
	}

	clock_t begin = clock();

	const int NFITPARAMETERS = 19;
	const double THRESHOLDBACKGROUND = 0.7;
	const int NDIMENSIONS = 4;
	const int NPOOLS = 7;

	const char *matlabFileIn = NULL;
	const char *cestDataIn = NULL;
	const char *frequencyOffsetIn = NULL;
	const char *frequencyOffsetOut = NULL;
	const char *poolSwitchOnOff = NULL;
	const char *matlabFileOut = NULL;
	char errorMessage[100];
	int Niterations;
	char parameterOption;
	ExtraData extraData;
	double arr[NPOOLS];
	extraData.singlePoolOn = arr;

	MATFile *matFileIn, *matFileOut;

	mxArray *pCestDataIn, *pFrequencyOffsetIn, *pFrequencyOffsetOut,
			*pFitParameters, *pPoolSwitchOnOff, *pMtPoolFit, *pNoe1PoolFit,
			*pNoe2PoolFit, *pCrPoolFit, *pAmidePoolFit, *pWaterPoolFit;

	double *cestDataInArray, *frequencyOffsetInArray, *frequencyOffsetOutArray,
			*poolSwitchOnOffArray;

	mwSize NfrequencyOffsetIn, NfrequencyOffsetOut, NpoolSwitchOnOff;

	const mwSize *dimensions;

	while ((parameterOption = getopt(argc, argv, "i:d:k:o:l:s:r:")) != EOF)
		switch (parameterOption) {
		case 'i':
			matlabFileIn = optarg;
			break;
		case 'd':
			cestDataIn = optarg;
			break;
		case 'k':
			frequencyOffsetIn = optarg;
			break;
		case 'o':
			matlabFileOut = optarg;
			break;
		case 'l':
			frequencyOffsetOut = optarg;
			break;
		case 's':
			poolSwitchOnOff = optarg;
			break;
		case 'r':
			Niterations = atoi(optarg);
			break;
		default:
			sprintf(errorMessage, "Unknown option: %s", optarg);
			error(errorMessage);
		}

	argc -= optind;
	argv += optind;

	char logFileName[100];
	strcpy(logFileName, matlabFileOut);
	strcat(logFileName, "_logfile.txt");

	FILE *logFile = fopen(logFileName, "w");
	fprintf(logFile, "%s\n", timeStamp());
	for (int count = 0; count < argc; count++) {
		fprintf(logFile, " %s ", argv[count]);
	}
	fclose(logFile);

// Open .mat file
	matFileIn = matOpen(matlabFileIn, "r");
	if (!matFileIn) {
		sprintf(errorMessage, "Cannot open %s", matlabFileIn);
		error(errorMessage);
	}

// Read CEST data
	pCestDataIn = matGetVariable(matFileIn, cestDataIn);
	cestDataInArray = mxGetDoubles(pCestDataIn);
	if (!pCestDataIn) {
		sprintf(errorMessage, "Cannot get %s", cestDataIn);
		error(errorMessage);
	}

// Read input offsets
	pFrequencyOffsetIn = matGetVariable(matFileIn, frequencyOffsetIn);
	NfrequencyOffsetIn = mxGetNumberOfElements(pFrequencyOffsetIn);
	frequencyOffsetInArray = mxGetDoubles(pFrequencyOffsetIn);
	if (!pFrequencyOffsetIn) {
		sprintf(errorMessage, "Cannot get %s", frequencyOffsetIn);
		error(errorMessage);
	}

// Read output offsets
	pFrequencyOffsetOut = matGetVariable(matFileIn, frequencyOffsetOut);
	NfrequencyOffsetOut = mxGetNumberOfElements(pFrequencyOffsetOut);
	frequencyOffsetOutArray = mxGetDoubles(pFrequencyOffsetOut);
	if (!pFrequencyOffsetOut) {
		sprintf(errorMessage, "Cannot get %s", frequencyOffsetOut);
		error(errorMessage);
	}

// Read pool switch
	pPoolSwitchOnOff = matGetVariable(matFileIn, poolSwitchOnOff);
	NpoolSwitchOnOff = mxGetNumberOfElements(pPoolSwitchOnOff);
	poolSwitchOnOffArray = mxGetDoubles(pPoolSwitchOnOff);
	if (!pPoolSwitchOnOff) {
		sprintf(errorMessage, "Cannot get %s", poolSwitchOnOff);
		error(errorMessage);
	}

	extraData.NoffsetsIn = NfrequencyOffsetIn;
	extraData.NoffsetsOut = NfrequencyOffsetOut;
	extraData.offsetsIn = frequencyOffsetInArray;
	extraData.offsetsOut = frequencyOffsetOutArray;
	extraData.poolSwitch = poolSwitchOnOffArray;
	dimensions = mxGetDimensions(pCestDataIn);

// Dimensions of CEST data
	int dim1 = dimensions[0];
	int dim2 = dimensions[1];
	int dim3 = dimensions[2];
	int dim4 = dimensions[3];

	int statusFitParameters, statusMtPoolFit, statusNoe1PoolFit,
			statusNoe2PoolFit, statusCrPoolFit, statusAmidePoolFit,
			statusWaterPoolFit;

// Dynamic array allocation for input CEST spectra
	double (*dataInArray)[dim2][dim3][dim4];
	dataInArray = malloc(dim1 * sizeof *dataInArray);

// Dynamic array allocation for fitting parameters
	double (*fitParametersArray)[dim2][dim3][NFITPARAMETERS];
	fitParametersArray = malloc(dim1 * sizeof *fitParametersArray);

// Dynamic array allocation for fitted MT pool spectra
	double (*fitMtPoolArray)[dim2][dim3][NfrequencyOffsetOut];
	fitMtPoolArray = malloc(dim1 * sizeof *fitMtPoolArray);

// Dynamic array allocation for fitted NOE1 pool spectra
	double (*fitNoe1PoolArray)[dim2][dim3][NfrequencyOffsetOut];
	fitNoe1PoolArray = malloc(dim1 * sizeof *fitNoe1PoolArray);

// Dynamic array allocation for fitted NOE2 pool spectra
	double (*fitNoe2PoolArray)[dim2][dim3][NfrequencyOffsetOut];
	fitNoe2PoolArray = malloc(dim1 * sizeof *fitNoe2PoolArray);

// Dynamic array allocation for fitted Cr pool spectra
	double (*fitCrPoolArray)[dim2][dim3][NfrequencyOffsetOut];
	fitCrPoolArray = malloc(dim1 * sizeof *fitCrPoolArray);

// Dynamic array allocation for fitted Amide pool spectra
	double (*fitAmidePoolArray)[dim2][dim3][NfrequencyOffsetOut];
	fitAmidePoolArray = malloc(dim1 * sizeof *fitAmidePoolArray);

// Dynamic array allocation for fitted Water pool spectra
	double (*fitWaterPoolArray)[dim2][dim3][NfrequencyOffsetOut];
	fitWaterPoolArray = malloc(dim1 * sizeof *fitWaterPoolArray);

// Creating output .mat file
	matFileOut = matOpen(matlabFileOut, "w");

	if (!matFileOut) {
		sprintf(errorMessage, "Error creating file %s", matlabFileOut);
		error(errorMessage);
	}

// MATLAB array for fitting parameters
	const mwSize dimsFitParameters[NDIMENSIONS] = { dim1, dim2, dim3,
			NFITPARAMETERS };
	pFitParameters = mxCreateNumericArray(NDIMENSIONS, dimsFitParameters,
			mxDOUBLE_CLASS, mxREAL);

	if (!pFitParameters) {
		sprintf(errorMessage,
				"Unable to create mxArray. %s : Out of memory on line %d",
				__FILE__, __LINE__);
		error(errorMessage);
	}

// MATLAB array for fitted MT pool spectra
	const mwSize dimsMtPoolFit[NDIMENSIONS] = { dim1, dim2, dim3,
			NfrequencyOffsetOut };
	pMtPoolFit = mxCreateNumericArray(NDIMENSIONS, dimsMtPoolFit,
			mxDOUBLE_CLASS, mxREAL);

	if (!pMtPoolFit) {
		sprintf(errorMessage,
				"Unable to create mxArray. %s : Out of memory on line %d",
				__FILE__, __LINE__);
		error(errorMessage);
	}

// MATLAB array for fitted NOE1 pool spectra
	const mwSize dimsNoe1PoolFit[NDIMENSIONS] = { dim1, dim2, dim3,
			NfrequencyOffsetOut };
	pNoe1PoolFit = mxCreateNumericArray(NDIMENSIONS, dimsNoe1PoolFit,
			mxDOUBLE_CLASS, mxREAL);

	if (!pNoe1PoolFit) {
		sprintf(errorMessage,
				"Unable to create mxArray. %s : Out of memory on line %d",
				__FILE__, __LINE__);
		error(errorMessage);
	}

// MATLAB array for fitted NOE2 pool spectra
	const mwSize dimsNoe2PoolFit[NDIMENSIONS] = { dim1, dim2, dim3,
			NfrequencyOffsetOut };
	pNoe2PoolFit = mxCreateNumericArray(NDIMENSIONS, dimsNoe2PoolFit,
			mxDOUBLE_CLASS, mxREAL);

	if (!pNoe2PoolFit) {
		sprintf(errorMessage,
				"Unable to create mxArray. %s : Out of memory on line %d",
				__FILE__, __LINE__);
		error(errorMessage);
	}

// MATLAB array for fitted Cr pool spectra
	const mwSize dimsCrPoolFit[NDIMENSIONS] = { dim1, dim2, dim3,
			NfrequencyOffsetOut };
	pCrPoolFit = mxCreateNumericArray(NDIMENSIONS, dimsCrPoolFit,
			mxDOUBLE_CLASS, mxREAL);

	if (!pCrPoolFit) {
		sprintf(errorMessage,
				"Unable to create mxArray. %s : Out of memory on line %d",
				__FILE__, __LINE__);
		error(errorMessage);
	}

// MATLAB array for fitted Amide pool spectra
	const mwSize dimsAmidePoolFit[NDIMENSIONS] = { dim1, dim2, dim3,
			NfrequencyOffsetOut };
	pAmidePoolFit = mxCreateNumericArray(NDIMENSIONS, dimsAmidePoolFit,
			mxDOUBLE_CLASS, mxREAL);

	if (!pAmidePoolFit) {
		sprintf(errorMessage,
				"Unable to create mxArray. %s : Out of memory on line %d",
				__FILE__, __LINE__);
		error(errorMessage);
	}

// MATLAB array for fitted Water pool spectra
	const mwSize dimsWaterPoolFit[NDIMENSIONS] = { dim1, dim2, dim3,
			NfrequencyOffsetOut };
	pWaterPoolFit = mxCreateNumericArray(NDIMENSIONS, dimsWaterPoolFit,
			mxDOUBLE_CLASS, mxREAL);

	if (!pWaterPoolFit) {
		sprintf(errorMessage,
				"Unable to create mxArray. %s : Out of memory on line %d",
				__FILE__, __LINE__);
		error(errorMessage);
	}

// MATLAB 4d array
	size_t idx1 = 0;
	for (size_t l = 0; l < dim4; l++) {
		for (size_t k = 0; k < dim3; k++) {
			for (size_t j = 0; j < dim2; j++) {
				for (size_t i = 0; i < dim1; i++) {
					dataInArray[i][j][k][l] = cestDataInArray[idx1];
					idx1++;
				}
			}
		}
	}

	double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
	opts[0] = 1E-03;
	opts[1] = 1E-08;
	opts[2] = 1E-08;
	opts[3] = 1E-08;
	opts[4] = 1E-06;

	for (size_t l = 0; l < dim1; l++) {

		double progress = (l * 100.0) / dim1;
		fprintf(stderr, "**%.2f%%**", progress);

		for (size_t k = 0; k < dim2; k++) {

			for (size_t j = 0; j < dim3; j++) {

				double sum = 0;
				double spectrum[NfrequencyOffsetIn];

				for (size_t i = 0; i < dim4; i++) {
					spectrum[i] = 1.0 - dataInArray[l][k][j][i];
					sum = sum + spectrum[i];
				}

				if (sum > NfrequencyOffsetIn * THRESHOLDBACKGROUND) {
					// Don't waste time fitting "bad" data
					continue;
				}

				double p[NFITPARAMETERS];
				p[0] = 0.2; // constant

				p[1] = 0.1; // MT amp
				p[2] = 15.0; // MT fwhm
				p[3] = -2.5; // MT cs

				p[4] = 0.2; // NOE1 amp
				p[5] = 3.0; // NOE1 fwhm
				p[6] = -3.5; // NOE1 cs

				p[7] = 0.02; // NOE2 amp
				p[8] = 3.0; // NOE2 fwhm
				p[9] = -1.6; // NOE2 cs

				p[10] = 0.05; // Cr amp
				p[11] = 1.0; // Cr fwhm
				p[12] = 2.0; // Cr cs

				p[13] = 0.025; // Amide amp
				p[14] = 1.0; // Amide fwhm
				p[15] = 3.5; // Amide cs

				p[16] = 0.9; // Water amp
				p[17] = 1.0; // Water fwhm
				p[18] = 0.0; // Water cs

				double lb[NFITPARAMETERS];
				lb[0] = 0.0; // constant

				lb[1] = 0.0; // MT amp
				lb[2] = 10.0; // MT fwhm
				lb[3] = -3.5; // MT cs

				lb[4] = 0.0; // NOE1 amp
				lb[5] = 0.5; // NOE1 fwhm
				lb[6] = -4.0; // NOE1 cs

				lb[7] = 0.0; // NOE2 amp
				lb[8] = 0.5; // NOE2 fwhm
				lb[9] = -1.8; // NOE2 cs

				lb[10] = 0.0; // Cr amp
				lb[11] = 0.3; // Cr fwhm
				lb[12] = 1.7; // Cr cs

				lb[13] = 0.0; // Amide amp
				lb[14] = 0.3; // Amide fwhm
				lb[15] = 3.3; // Amide cs

				lb[16] = 0.2; // DS amp
				lb[17] = 0.2; // DS fwhm
				lb[18] = -1.0; // DS cs

				double ub[NFITPARAMETERS];
				ub[0] = 0.5; // constant

				ub[1] = 1.0; // MT amp
				ub[2] = 100; // MT fwhm
				ub[3] = -2.0; // MT cs

				ub[4] = 0.8; // NOE1 amp
				ub[5] = 5.0; // NOE1 fwhm
				ub[6] = -3.2; // NOE1 cs

				ub[7] = 0.25; // NOE2 amp
				ub[8] = 3.0; // NOE2 fwhm
				ub[9] = -1.4; // NOE2 cs

				ub[10] = 0.3; // Cr amp
				ub[11] = 1.0; // Cr fwhm
				ub[12] = 2.2; // Cr cs

				ub[13] = 0.3; // Amide amp
				ub[14] = 2.5; // Amide fwhm
				ub[15] = 3.7; // Amide cs

				ub[16] = 1.0; // Water amp
				ub[17] = 3.0; // Water fwhm
				ub[18] = 1.0; // Water cs

				// Fitting algorithm NOT thread safe
				int ret = dlevmar_bc_dif(findFullModelFitParameters, p,
						spectrum, NFITPARAMETERS, NfrequencyOffsetIn, lb, ub,
						NULL, Niterations, opts, info, NULL, NULL,
						(void*) &extraData); // no Jacobian

				double singlePoolFittedSpectrum[NPOOLS][NfrequencyOffsetOut];
				for (int poolNo = 0; poolNo < NPOOLS; poolNo++) {
					memset(extraData.singlePoolOn, 0, NPOOLS * sizeof(double));
					(extraData.singlePoolOn)[poolNo] = 1;
					calculateArbitraryPoolFit(p,
							singlePoolFittedSpectrum[poolNo],
							(void*) &extraData);
				}

				for (size_t parNo = 0; parNo < NFITPARAMETERS; parNo++) {
					fitParametersArray[l][k][j][parNo] = p[parNo];
				}

				for (size_t offsetNo = 0; offsetNo < NfrequencyOffsetOut;
						offsetNo++) {
					fitMtPoolArray[l][k][j][offsetNo] =
							singlePoolFittedSpectrum[MT][offsetNo];
					fitNoe1PoolArray[l][k][j][offsetNo] =
							singlePoolFittedSpectrum[NOE1][offsetNo];
					fitNoe2PoolArray[l][k][j][offsetNo] =
							singlePoolFittedSpectrum[NOE2][offsetNo];
					fitCrPoolArray[l][k][j][offsetNo] =
							singlePoolFittedSpectrum[CR][offsetNo];
					fitAmidePoolArray[l][k][j][offsetNo] =
							singlePoolFittedSpectrum[AMIDE][offsetNo];
					fitWaterPoolArray[l][k][j][offsetNo] =
							singlePoolFittedSpectrum[WATER][offsetNo];
				}
			}
		}
	}

// Copy fitting pars into MATLAB array
	double *startFitParameters = (double*) mxGetDoubles(pFitParameters);
	size_t bytesFitParameters = (dim1 * dim2 * dim3 * NFITPARAMETERS)
			* mxGetElementSize(pFitParameters);
	memcpy(startFitParameters, fitParametersArray, bytesFitParameters);

// Put fitting pars array into .mat file
	statusFitParameters = matPutVariable(matFileOut, "fitPars", pFitParameters);

	if (statusFitParameters) {
		sprintf(errorMessage, "%s :  Error using matPutVariable on line %d\n",
		__FILE__, __LINE__);
		error(errorMessage);
	}

// Copy fitted MT pool spectra into MATLAB array
	double *startMtPoolFit = (double*) mxGetDoubles(pMtPoolFit);
	size_t bytesMtPoolFit = (dim1 * dim2 * dim3 * NfrequencyOffsetOut)
			* mxGetElementSize(pMtPoolFit);
	memcpy(startMtPoolFit, fitMtPoolArray, bytesMtPoolFit);

// Put fitted MT pool spectra into .mat file
	statusMtPoolFit = matPutVariable(matFileOut, "fitMtPool", pMtPoolFit);

	if (statusMtPoolFit) {
		sprintf(errorMessage, "%s :  Error using matPutVariable on line %d\n",
		__FILE__, __LINE__);
		error(errorMessage);
	}

// Copy fitted NOE1 pool spectra into MATLAB array
	double *startNoe1PoolFit = (double*) mxGetDoubles(pNoe1PoolFit);
	size_t bytesNoe1PoolFit = (dim1 * dim2 * dim3 * NfrequencyOffsetOut)
			* mxGetElementSize(pNoe1PoolFit);
	memcpy(startNoe1PoolFit, fitNoe1PoolArray, bytesNoe1PoolFit);

// Put fitted NOE1 pool spectra into .mat file
	statusNoe1PoolFit = matPutVariable(matFileOut, "fitNoe1Pool", pNoe1PoolFit);

	if (statusNoe1PoolFit) {
		sprintf(errorMessage, "%s :  Error using matPutVariable on line %d\n",
		__FILE__, __LINE__);
		error(errorMessage);
	}

// Copy fitted NOE2 pool spectra into MATLAB array
	double *startNoe2PoolFit = (double*) mxGetDoubles(pNoe2PoolFit);
	size_t bytesNoe2PoolFit = (dim1 * dim2 * dim3 * NfrequencyOffsetOut)
			* mxGetElementSize(pNoe2PoolFit);
	memcpy(startNoe2PoolFit, fitNoe2PoolArray, bytesNoe2PoolFit);

// Put fitted NOE2 pool spectra into .mat file
	statusNoe2PoolFit = matPutVariable(matFileOut, "fitNoe2Pool", pNoe2PoolFit);

	if (statusNoe2PoolFit) {
		sprintf(errorMessage, "%s :  Error using matPutVariable on line %d\n",
		__FILE__, __LINE__);
		error(errorMessage);
	}

// Copy fitted Cr pool spectra into MATLAB array
	double *startCrPoolFit = (double*) mxGetDoubles(pCrPoolFit);
	size_t bytesCrPoolFit = (dim1 * dim2 * dim3 * NfrequencyOffsetOut)
			* mxGetElementSize(pCrPoolFit);
	memcpy(startCrPoolFit, fitCrPoolArray, bytesCrPoolFit);

// Put fitted Crpool spectra into .mat file
	statusCrPoolFit = matPutVariable(matFileOut, "fitCrPool", pCrPoolFit);

	if (statusCrPoolFit) {
		sprintf(errorMessage, "%s :  Error using matPutVariable on line %d\n",
		__FILE__, __LINE__);
		error(errorMessage);
	}

// Copy fitted Amide pool spectra into MATLAB array
	double *startAmidePoolFit = (double*) mxGetDoubles(pAmidePoolFit);
	size_t bytesAmidePoolFit = (dim1 * dim2 * dim3 * NfrequencyOffsetOut)
			* mxGetElementSize(pAmidePoolFit);
	memcpy(startAmidePoolFit, fitAmidePoolArray, bytesAmidePoolFit);

// Put fitted Amide pool spectra into .mat file
	statusAmidePoolFit = matPutVariable(matFileOut, "fitAmidePool",
			pAmidePoolFit);

	if (statusAmidePoolFit) {
		sprintf(errorMessage, "%s :  Error using matPutVariable on line %d\n",
		__FILE__, __LINE__);
		error(errorMessage);
	}

// Copy fitted Water pool spectra into MATLAB array
	double *startWaterPoolFit = (double*) mxGetDoubles(pWaterPoolFit);
	size_t bytesWaterPoolFit = (dim1 * dim2 * dim3 * NfrequencyOffsetOut)
			* mxGetElementSize(pWaterPoolFit);
	memcpy(startWaterPoolFit, fitWaterPoolArray, bytesWaterPoolFit);

// Put fitted Water pool spectra into .mat file
	statusWaterPoolFit = matPutVariable(matFileOut, "fitWaterPool",
			pWaterPoolFit);

	if (statusWaterPoolFit) {
		sprintf(errorMessage, "%s :  Error using matPutVariable on line %d\n",
		__FILE__, __LINE__);
		error(errorMessage);
	}

// Clean up
	mxDestroyArray(pCestDataIn);
	mxDestroyArray(pFitParameters);
	mxDestroyArray(pMtPoolFit);
	mxDestroyArray(pNoe1PoolFit);
	mxDestroyArray(pNoe2PoolFit);
	mxDestroyArray(pCrPoolFit);
	mxDestroyArray(pAmidePoolFit);
	mxDestroyArray(pWaterPoolFit);

	free(dataInArray);
	free(fitParametersArray);
	free(fitMtPoolArray);
	free(fitNoe1PoolArray);
	free(fitNoe2PoolArray);
	free(fitCrPoolArray);
	free(fitAmidePoolArray);
	free(fitWaterPoolArray);

	if (matClose(matFileIn)) {
		sprintf(errorMessage, "Error closing file %s\n", matlabFileIn);
		error(errorMessage);
	}

// Stop timer
	clock_t end = clock();
	double cpu_time_used = (double) (end - begin) / CLOCKS_PER_SEC;
	printf("\nFitting done in [s] %.10g \n: ", cpu_time_used);

	return 0;
}



