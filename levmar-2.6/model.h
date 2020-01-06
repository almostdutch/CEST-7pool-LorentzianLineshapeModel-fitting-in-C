
// MT - magnetization transfer pool at -2.4 ppm
void calculateMtPoolFit(double *p, double *result, void *data);

// NOE1 - nuclear overhauser enhancement at -3.5 ppm
void calculateNoe1PoolFit(double *p, double *result, void *data);

// NOE2 - nuclear overhauser enhancement at -1.6 ppm
void calculateNoe2poolFit(double *p, double *result, void *data);

// Cr - creatine (and other amine CEST effects) 2.0 ppm
void calculateCrPoolFit(double *p, double *result, void *data);

// Amide - (poly)-peptides at 3.5 ppm
void calculateAmidePoolFit(double *p, double *result, void *data);

// Direct water saturation at 0 ppm
void calculateWaterPoolFit(double *p, double *result, void *data);

// Cost function for data fitting
void findFullModelFitParameters(double *p, double *x, int m, int n, void* data);

