#define MAX(a, b) ((a) > (b)  ? (a) : (b))
#define MIN(a, b) ((a) < (b)  ? (a) : (b))
#define EPS pow(10, -9)

void probSxSyCond(int *vx, int *vy, double *logpxy, double *logj, double *factj,
		  int nx, int ny, int nux, int nuy, int nuxy, int *ixy, int *iyx,
		  double combx, double comby, double *sprob, int *mmax, int nm);
double probSxSy(double *rm, double *sprob, int nm, int mmax);
double probSxSyEqr(double *rm, double *sprob, int nm, int mmax, double *cbinom);
int equalArr(int *arr1, int *arr2, int narr);
void vMax(int *base, int nbase, int sumlim, int *vmax);
