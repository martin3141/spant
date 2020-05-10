#include <math.h>
// #include <R_ext/Print.h>

// Find the index of each x[i] in a vector v sorted in ascending order,
// i.e., for i = 1, ..., m, returns
// 
//           min {0 <= j <= n: x[i] <= v[j]}
//
// where v[0] = -Inf.

void indx(double *x, int *m, double *v, int *n, int *ind)
{
  int left, right, mid;
  int i, j;

  for(i = 0; i < *m; i++) {
    if(x[i] < v[0]) {
      ind[i] = 0;
      continue;
    }
    if(x[i] >= v[*n-1]) {
      ind[i] = *n;
      continue;
    }
    left = 0; 
    right = *n - 1;            // always x[i] < x[right]
    while(1) {
      if(left >= right - 1) {
	ind[i] = left + 1;
	break;
      }
      mid = rint((left + right) * 0.5);
      if(x[i] >= v[mid]) {     // to find the rightmost tied value
	left = mid;
      }
      else {
	right = mid;
      }
    }
  }
}

// Returns either row or column maxima

void matMaxs(double *x, int *n, int *m, double *v, int *dim)
{
  int i, j;
  
  if(*dim == 1) {         // rows
    for(i = 0; i < *n; i++) {
      v[i] = x[i];
      if(*m > 0) {
	for(j = 1; j < *m; j++) {
	  if(x[j*(*n)+i] > v[i]) {
	    v[i] = x[j*(*n)+i];
	  }
	}
      }
    }
  }
  else {                  // columns
    for(j = 0; j < *m; j++) {
      v[j] = x[j*(*n)];
      if(*n > 0) {
	for(i = 1; i < *n; i++) {
	  if(x[j*(*n)+i] > v[j]) {
	    v[j] = x[j*(*n)+i];
	  }
	}
      }
    }
  }
}
