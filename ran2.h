#ifndef RAN2_H
#define RAN2_H
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-15
#define RNMX (1.0-EPS)


/*Long period (>2*10^18) random number generator of L'Ecuyer 
 * with Bays-Durham shuffle and added safeguards. 
 * Ruturns a uniform random deviate between 0.0 and 1.0(exclusive of the endpoint values).
 * Call with idum a negative integer to initialize; thereafter,
 * do not alter idum between succesive deviates in a sequence.
 * RNMX should approximate the largest floating value that is less than 1.*/
float ran2(long *idum);

#endif
