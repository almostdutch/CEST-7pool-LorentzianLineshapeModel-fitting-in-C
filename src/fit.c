#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <errno.h>
#include <time.h>
#include <math.h>
#include "multiLorentzian.h"

void error(char *msg) {
	fprintf(stderr, "%s: %s\n", msg, strerror(errno));
	exit(1);
}

char* timeStamp() {
	time_t t;
	time(&t);
	return asctime(localtime(&t));
}

int catchSignal(int sig, void (*handler)(int)) {
	struct sigaction action;
	action.sa_handler = handler;
	sigemptyset(&action.sa_mask);
	action.sa_flags = 0;
	return sigaction(sig, &action, NULL);
}

void interrupt(int sig) {
	puts("Terminated with Ctrl + C");
	puts("Exiting the program");
	exit(1);
}

void findFullModelFitParameters(double *p, double *spectrum, int NFITPARAMETERS,
		int NfrequencyOffsetIn, void *data) {
	ExtraData *dptr;
	dptr = (ExtraData*) data;

	double singleOffset;
	for (int i = 0; i < NfrequencyOffsetIn; ++i) {
		singleOffset = dptr->offsetsIn[i];

		spectrum[i] = p[0] * dptr->poolSwitch[0]
				+ p[1] * (1 / (1 + pow((singleOffset - p[3]) / (p[2]), 2)))
						* dptr->poolSwitch[1]
				+ p[4] * (1 / (1 + pow((singleOffset - p[6]) / (p[5]), 2)))
						* dptr->poolSwitch[2]
				+ p[7] * (1 / (1 + pow((singleOffset - p[9]) / (p[8]), 2)))
						* dptr->poolSwitch[3]
				+ p[10] * (1 / (1 + pow((singleOffset - p[12]) / (p[11]), 2)))
						* dptr->poolSwitch[4]
				+ p[13] * (1 / (1 + pow((singleOffset - p[15]) / (p[14]), 2)))
						* dptr->poolSwitch[5]
				+ p[16] * (1 / (1 + pow((singleOffset - p[18]) / (p[17]), 2)))
						* dptr->poolSwitch[6];
	}
}

void calculateArbitraryPoolFit(double *p, double *spectrum, void *data) {
	ExtraData *dptr;
	dptr = (ExtraData*) data;

	double singleOffset;
	for (int i = 0; i < dptr->NoffsetsOut; ++i) {
		singleOffset = dptr->offsetsOut[i];

		spectrum[i] = p[0] * dptr->poolSwitch[0] * dptr->singlePoolOn[0]
				+ p[1] * (1 / (1 + pow((singleOffset - p[3]) / (p[2]), 2)))
						* dptr->poolSwitch[1] * dptr->singlePoolOn[1]
				+ p[4] * (1 / (1 + pow((singleOffset - p[6]) / (p[5]), 2)))
						* dptr->poolSwitch[2] * dptr->singlePoolOn[2]
				+ p[7] * (1 / (1 + pow((singleOffset - p[9]) / (p[8]), 2)))
						* dptr->poolSwitch[3] * dptr->singlePoolOn[3]
				+ p[10] * (1 / (1 + pow((singleOffset - p[12]) / (p[11]), 2)))
						* dptr->poolSwitch[4] * dptr->singlePoolOn[4]
				+ p[13] * (1 / (1 + pow((singleOffset - p[15]) / (p[14]), 2)))
						* dptr->poolSwitch[5] * dptr->singlePoolOn[5]
				+ p[16] * (1 / (1 + pow((singleOffset - p[18]) / (p[17]), 2)))
						* dptr->poolSwitch[6] * dptr->singlePoolOn[6];
	}
}
