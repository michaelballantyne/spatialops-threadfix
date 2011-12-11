#include <stdio.h>
#include <stdlib.h>
#define MAX_THREADS_PER_BLOCK 256

__global__ void _vecAdd(float* vec1, float* vec2, float* returnVec,
		int vecSize) {
	int i;

	i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < vecSize) {
		returnVec[i] = vec1[i] + vec2[i];
	}
}

__global__ void _vecMul(float* vec1, float* vec2, float* returnVec,
		int vecSize) {
	int i;

	i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < vecSize) {
		returnVec[i] = vec1[i] * vec2[i];
	}
}

//Device vector add implementation
__host__ float* vecAdd(float* vec1, float* vec2, int vecSize) {
	float *d_vec1, *d_vec2;
	float *d_returnVec;
	float *h_returnVec;
	int vecMemSize;
	int blockWidth;
	int numOfBlocks;

	blockWidth =
			MAX_THREADS_PER_BLOCK < vecSize ? MAX_THREADS_PER_BLOCK : vecSize;
	numOfBlocks = vecSize / blockWidth + (vecSize % blockWidth == 0 ? 0 : 1);

	dim3 dimBlock(blockWidth, 1);
	dim3 dimGrid(numOfBlocks, 1);

	vecMemSize = vecSize * sizeof(float);
	h_returnVec = (float*) malloc(vecMemSize);
	cudaMalloc((void**) &d_vec1, vecMemSize);
	cudaMalloc((void**) &d_vec2, vecMemSize);
	cudaMalloc((void**) &d_returnVec, vecMemSize);

	cudaMemcpy(d_vec1, vec1, vecMemSize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_vec2, vec2, vecMemSize, cudaMemcpyHostToDevice);

	_vecAdd<<<dimGrid, dimBlock, 0, 0>>>(d_vec1, d_vec2, d_returnVec, vecSize);
	cudaThreadSynchronize();

	cudaError err;
	err = cudaGetLastError();
	if (cudaSuccess != err) {
		fprintf(stderr, "Function call failed!\n");
		fprintf(stderr, "(Why?) : %s\n", cudaGetErrorString(err));
	}

	cudaMemcpy(h_returnVec, d_returnVec, vecMemSize, cudaMemcpyDeviceToHost);
	cudaFree(d_vec1);
	cudaFree(d_vec2);
	cudaFree(d_returnVec);

	return h_returnVec;
}

//Device vector mul implementation
__host__ float* vecMul(float* vec1, float* vec2, int vecSize) {
	float *d_vec1, *d_vec2;
	float *d_returnVec;
	float *h_returnVec;
	int vecMemSize;
	int blockWidth;
	int numOfBlocks;

	blockWidth =
			MAX_THREADS_PER_BLOCK < vecSize ? MAX_THREADS_PER_BLOCK : vecSize;
	numOfBlocks = vecSize / blockWidth + (vecSize % blockWidth == 0 ? 0 : 1);

	dim3 dimBlock(blockWidth, 1);
	dim3 dimGrid(numOfBlocks, 1);

	vecMemSize = vecSize * sizeof(float);
	h_returnVec = (float*) malloc(vecMemSize);
	cudaMalloc((void**) &d_vec1, vecMemSize);
	cudaMalloc((void**) &d_vec2, vecMemSize);
	cudaMalloc((void**) &d_returnVec, vecMemSize);

	cudaMemcpy(d_vec1, vec1, vecMemSize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_vec2, vec2, vecMemSize, cudaMemcpyHostToDevice);

	_vecMul<<<dimGrid, dimBlock, 0, 0>>>(d_vec1, d_vec2, d_returnVec, vecSize);
	cudaThreadSynchronize();

	cudaMemcpy(h_returnVec, d_returnVec, vecMemSize, cudaMemcpyDeviceToHost);
	cudaFree(d_vec1);
	cudaFree(d_vec2);
	cudaFree(d_returnVec);

	return h_returnVec;
}

int func_add(float *x, float *y, int sz) {
	int i;
	float *a;
	a = (float *) malloc(sizeof(float) * sz);
	if (!a) {
		printf("memory allocation error\n");
		exit(-1);
	}
	memcpy(a, x, sz * (sizeof(float)));

	x = vecAdd(x, y, sz);

	for (i = 0; i < sz; i++) {
		if (x[i] != a[i] + y[i]) {
			return 0;
		}
	}

	free(a);
	return 1;
}

int func_mul(float *x, float *y, int sz) {
	int i;
	float *a;

	a = (float *) malloc(sizeof(float) * sz);
	if (!a) {
		printf("memory allocation error\n");
		exit(-1);
	}
	memcpy(a, x, sz * (sizeof(float)));

	x = vecMul(x, y, sz);

	for (i = 0; i < sz; i++) {
		if (x[i] != a[i] * y[i]) {
			return 0;
		}
	}

	free(a);
	return 1;
}
