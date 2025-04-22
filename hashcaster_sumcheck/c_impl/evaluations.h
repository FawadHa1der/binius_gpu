#ifndef EVALUATION_H
#define EVALUATION_H

#include "mle_poly.h"

typedef Points Evaluations;

void twist_evals(Evaluations *evals);
void untwist_evals(Evaluations *evals);

#endif