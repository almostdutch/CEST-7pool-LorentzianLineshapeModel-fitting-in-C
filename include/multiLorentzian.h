#ifndef _FIT_H_
#define _FIT_H_

enum Pools {
	CONSTANT, MT, NOE1, NOE2, CR, AMIDE, WATER
};

typedef struct {
	double *offsetsIn;
	int NoffsetsIn;
	double *offsetsOut;
	int NoffsetsOut;
	double *poolSwitch;
	double *singlePoolOn;
} ExtraData;

void error(char *);
char* timeStamp(void);
int catchSignal(int, void (*handler)(int));
void interrupt(int);
void findFullModelFitParameters(double *, double *, int, int, void *);
void calculateArbitraryPoolFit(double *, double *, void *);

#endif
