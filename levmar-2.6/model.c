#include "model.h"

void calculateMtPoolFit(double *p, double *result, void *data) {
	struct pass_offset *dptr;
	dptr = (struct pass_offset *) data;

	double x[dptr->n_offset_OUT];
	double w;
	for (int i = 0; i < dptr->n_offset_OUT; ++i) {
		w = dptr->offset_OUT[i];

		x[i] = 0 * p[0] * dptr->pool_switch[0]
				+ 1 * p[1] * (1 / (1 + pow((w - p[3]) / (p[2]), 2)))
						* dptr->pool_switch[1]
				+ 0 * p[4] * (1 / (1 + pow((w - p[6]) / (p[5]), 2)))
						* dptr->pool_switch[2]
				+ 0 * p[7] * (1 / (1 + pow((w - p[9]) / (p[8]), 2)))
						* dptr->pool_switch[3]
				+ 0 * p[10] * (1 / (1 + pow((w - p[12]) / (p[11]), 2)))
						* dptr->pool_switch[4]
				+ 0 * p[13] * (1 / (1 + pow((w - p[15]) / (p[14]), 2)))
						* dptr->pool_switch[5]
				+ 0 * p[16] * (1 / (1 + pow((w - p[18]) / (p[17]), 2)))
						* dptr->pool_switch[6];

		result[i] = x[i];
	}
	return;
}

void calculateNoe1PoolFit(double *p, double *result, void *data) {
	struct pass_offset *dptr;
	dptr = (struct pass_offset *) data;

	double x[dptr->n_offset_OUT];
	double w;
	for (int i = 0; i < dptr->n_offset_OUT; ++i) {
		w = dptr->offset_OUT[i];

		x[i] = 0 * p[0] * dptr->pool_switch[0]
				+ 0 * p[1] * (1 / (1 + pow((w - p[3]) / (p[2]), 2)))
						* dptr->pool_switch[1]
				+ 1 * p[4] * (1 / (1 + pow((w - p[6]) / (p[5]), 2)))
						* dptr->pool_switch[2]
				+ 0 * p[7] * (1 / (1 + pow((w - p[9]) / (p[8]), 2)))
						* dptr->pool_switch[3]
				+ 0 * p[10] * (1 / (1 + pow((w - p[12]) / (p[11]), 2)))
						* dptr->pool_switch[4]
				+ 0 * p[13] * (1 / (1 + pow((w - p[15]) / (p[14]), 2)))
						* dptr->pool_switch[5]
				+ 0 * p[16] * (1 / (1 + pow((w - p[18]) / (p[17]), 2)))
						* dptr->pool_switch[6];

		result[i] = x[i];
	}
	return;
}

void calculateNoe2poolFit(double *p, double *result, void *data) {
	struct pass_offset *dptr;
	dptr = (struct pass_offset *) data;

	double x[dptr->n_offset_OUT];
	double w;
	for (int i = 0; i < dptr->n_offset_OUT; ++i) {
		w = dptr->offset_OUT[i];

		x[i] = 0 * p[0] * dptr->pool_switch[0]
				+ 0 * p[1] * (1 / (1 + pow((w - p[3]) / (p[2]), 2)))
						* dptr->pool_switch[1]
				+ 0 * p[4] * (1 / (1 + pow((w - p[6]) / (p[5]), 2)))
						* dptr->pool_switch[2]
				+ 1 * p[7] * (1 / (1 + pow((w - p[9]) / (p[8]), 2)))
						* dptr->pool_switch[3]
				+ 0 * p[10] * (1 / (1 + pow((w - p[12]) / (p[11]), 2)))
						* dptr->pool_switch[4]
				+ 0 * p[13] * (1 / (1 + pow((w - p[15]) / (p[14]), 2)))
						* dptr->pool_switch[5]
				+ 0 * p[16] * (1 / (1 + pow((w - p[18]) / (p[17]), 2)))
						* dptr->pool_switch[6];

		result[i] = x[i];
	}
	return;
}

void calculateCrPoolFit(double *p, double *result, void *data) {
	struct pass_offset *dptr;
	dptr = (struct pass_offset *) data;

	double x[dptr->n_offset_OUT];
	double w;
	for (int i = 0; i < dptr->n_offset_OUT; ++i) {
		w = dptr->offset_OUT[i];

		x[i] = 0 * p[0] * dptr->pool_switch[0]
				+ 0 * p[1] * (1 / (1 + pow((w - p[3]) / (p[2]), 2)))
						* dptr->pool_switch[1]
				+ 0 * p[4] * (1 / (1 + pow((w - p[6]) / (p[5]), 2)))
						* dptr->pool_switch[2]
				+ 0 * p[7] * (1 / (1 + pow((w - p[9]) / (p[8]), 2)))
						* dptr->pool_switch[3]
				+ 1 * p[10] * (1 / (1 + pow((w - p[12]) / (p[11]), 2)))
						* dptr->pool_switch[4]
				+ 0 * p[13] * (1 / (1 + pow((w - p[15]) / (p[14]), 2)))
						* dptr->pool_switch[5]
				+ 0 * p[16] * (1 / (1 + pow((w - p[18]) / (p[17]), 2)))
						* dptr->pool_switch[6];

		result[i] = x[i];
	}
	return;
}

void calculateAmidePoolFit(double *p, double *result, void *data) {
	struct pass_offset *dptr;
	dptr = (struct pass_offset *) data;

	double x[dptr->n_offset_OUT];
	double w;
	for (int i = 0; i < dptr->n_offset_OUT; ++i) {
		w = dptr->offset_OUT[i];

		x[i] = 0 * p[0] * dptr->pool_switch[0]
				+ 0 * p[1] * (1 / (1 + pow((w - p[3]) / (p[2]), 2)))
						* dptr->pool_switch[1]
				+ 0 * p[4] * (1 / (1 + pow((w - p[6]) / (p[5]), 2)))
						* dptr->pool_switch[2]
				+ 0 * p[7] * (1 / (1 + pow((w - p[9]) / (p[8]), 2)))
						* dptr->pool_switch[3]
				+ 0 * p[10] * (1 / (1 + pow((w - p[12]) / (p[11]), 2)))
						* dptr->pool_switch[4]
				+ 1 * p[13] * (1 / (1 + pow((w - p[15]) / (p[14]), 2)))
						* dptr->pool_switch[5]
				+ 0 * p[16] * (1 / (1 + pow((w - p[18]) / (p[17]), 2)))
						* dptr->pool_switch[6];

		result[i] = x[i];
	}
	return;
}

void calculateWaterPoolFit(double *p, double *result, void *data) {
	struct pass_offset *dptr;
	dptr = (struct pass_offset *) data;

	double x[dptr->n_offset_OUT];
	double w;
	for (int i = 0; i < dptr->n_offset_OUT; ++i) {
		w = dptr->offset_OUT[i];

		x[i] = 0 * p[0] * dptr->pool_switch[0]
				+ 0 * p[1] * (1 / (1 + pow((w - p[3]) / (p[2]), 2)))
						* dptr->pool_switch[1]
				+ 0 * p[4] * (1 / (1 + pow((w - p[6]) / (p[5]), 2)))
						* dptr->pool_switch[2]
				+ 0 * p[7] * (1 / (1 + pow((w - p[9]) / (p[8]), 2)))
						* dptr->pool_switch[3]
				+ 0 * p[10] * (1 / (1 + pow((w - p[12]) / (p[11]), 2)))
						* dptr->pool_switch[4]
				+ 0 * p[13] * (1 / (1 + pow((w - p[15]) / (p[14]), 2)))
						* dptr->pool_switch[5]
				+ 1 * p[16] * (1 / (1 + pow((w - p[18]) / (p[17]), 2)))
						* dptr->pool_switch[6];

		result[i] = x[i];
	}
	return;
}

void findFullModelFitParameters(double *p, double *x, int m, int n, void *data) {
	struct pass_offset *dptr;
	dptr = (struct pass_offset *) data;

	double w;
	for (int i = 0; i < dptr->n_offset_IN; ++i) {
		w = dptr->offset_IN[i];

		x[i] = p[0] * dptr->pool_switch[0]
				+ p[1] * (1 / (1 + pow((w - p[3]) / (p[2]), 2)))
						* dptr->pool_switch[1]
				+ p[4] * (1 / (1 + pow((w - p[6]) / (p[5]), 2)))
						* dptr->pool_switch[2]
				+ p[7] * (1 / (1 + pow((w - p[9]) / (p[8]), 2)))
						* dptr->pool_switch[3]
				+ p[10] * (1 / (1 + pow((w - p[12]) / (p[11]), 2)))
						* dptr->pool_switch[4]
				+ p[13] * (1 / (1 + pow((w - p[15]) / (p[14]), 2)))
						* dptr->pool_switch[5]
				+ p[16] * (1 / (1 + pow((w - p[18]) / (p[17]), 2)))
						* dptr->pool_switch[6];
	}
	return;
}

