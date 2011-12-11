#ifndef __DVN_A1_IMP
#define __DVN_A1_IMP

//GPU vector addition wrapper
float* vecAdd( float* vec1, float* vec2, int vecSize );

	//GPU vector multiplication wrapper
float* vecMul(  float* vec1, float* vec2, int vecSize );

int func_add ( float *x, float *y, int sz);

int func_mul ( float *x, float *y, int sz);

#endif 
