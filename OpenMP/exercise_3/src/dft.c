#include "dft.h"

// DFT/IDFT routine
// idft: 1 direct DFT, -1 inverse IDFT (Inverse DFT)
int aDFT(int idft, double* xr, double* xi, double* Xr_o, double* Xi_o, int N)
{
	#pragma omp parallel
	{
		#pragma omp for 
		for (int k=0 ; k<N*N ; k++){
			// Real part of X[k]
			int n = k % N;
			Xr_o[k / N] += xr[n] * cos(n * (k / N) * PI2 / N) + idft*xi[n]*sin(n * (k / N) * PI2 / N);
			// Imaginary part of X[k]
			Xi_o[k / N] += -idft*xr[n] * sin(n * (k / N) * PI2 / N) + xi[n] * cos(n * (k / N) * PI2 / N);  
		} 
	}
		
	// normalize if you are doing IDFT
	if (idft == -1){
		for (int n = 0; n < N; n++) {
			Xr_o[n] /= N;
			Xi_o[n] /= N;
		}
	}
	return 1;
}

int DFT(int idft, double* xr, double* xi, double* Xr_o, double* Xi_o, int N)
{
	for (int k = 0 ; k < N; k++) {
		for (int n = 0; n < N; n++) {
			// Real part of X[k]
			Xr_o[k] += xr[n] * cos(n * k * PI2 / N) + idft * xi[n] * sin(n * k * PI2 / N);
			// Imaginary part of X[k]
			Xi_o[k] += -idft * xr[n] * sin(n * k * PI2 / N) + xi[n] * cos(n * k * PI2 / N);
		}
	}
		
	// normalize if you are doing IDFT
	if (idft == -1){
		for (int n = 0; n < N; n++) {
			Xr_o[n] /= N;
			Xi_o[n] /= N;
		}
	}
	return 1;
}