
#ifndef DISTSIZE
#error DISTSIZE has not been defined.
#endif

#define CONCATx(n,a,b) a ## n ## _ ## b
#define CONCAT(n,a,b) CONCATx(n,a,b)
#define dfun(b) CONCAT(DISTSIZE,dist,b)
#define cdfun(b) CONCAT(DISTSIZE,cond_dist,b)

