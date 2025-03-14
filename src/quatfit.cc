#include "sp_type.h"
/*
* S. K. Kearsley, "On the Orthogonal Transformation Used for
* Structural Comparisons", Acta Crystallographica Section A,
* 45, 208-210 (1989)
*/
/*=======================================================
	a - input: matrix to diagonalize
	v - output: eigenvectors
	d - output: eigenvalues
=========================================================*/

void jacobi(double a[][4], double *d, double v[][4]){
	const int nrot = 50;
	double onorm, dnorm;
	double b, dma, q, t, c, s;
	double atemp, vtemp, dtemp;
	int i, j, k, l;

	for(j=0; j<4; j++) {
		for(i=0; i<4; i++) v[i][j] = 0.0;
		v[j][j] = 1.0;
		d[j] = a[j][j];
	}

	for(l=0; l<nrot; l++){
		dnorm = 0.0;
		onorm = 0.0;
		for (j=0; j<4; j++) {
			dnorm = dnorm + fabs(d[j]);
			for (i=0; i<j; i++) {
				onorm = onorm + fabs(a[i][j]);
		  }
		}
		if((onorm / dnorm) <= 1.0e-12) break; //goto Exit_now;
//
		for (i=0; i<4; i++)
		for (j=i+1; j<4; j++){
			b = a[i][j];
			if(fabs(b) <= 0.0) break;
			dma = d[j] - d[i];
			if((fabs(dma) + fabs(b)) <=  fabs(dma)) t = b / dma;
			else {
				q = 0.5 * dma / b;
				t = 1.0 / (fabs(q) + sqrt(1.0+q*q));
				if(q < 0.0) t = -t;
			}
			c = 1.0/sqrt(t * t + 1.0);
			s = t * c;
			a[i][j] = 0.0;
			for (k = 0; k<i; k++) {
				atemp = c * a[k][i] - s * a[k][j];
				a[k][j] = s * a[k][i] + c * a[k][j];
				a[k][i] = atemp;
			}
			for (k=i+1; k<j; k++) {
				atemp = c * a[i][k] - s * a[k][j];
				a[k][j] = s * a[i][k] + c * a[k][j];
				a[i][k] = atemp;
				}
			for (k=j+1; k<4; k++) {
				atemp = c * a[i][k] - s * a[j][k];
				a[j][k] = s * a[i][k] + c * a[j][k];
				a[i][k] = atemp;
				}
			for (k=0; k<4; k++) {
				vtemp = c * v[k][i] - s * v[k][j];
				v[k][j] = s * v[k][i] + c * v[k][j];
				v[k][i] = vtemp;
			}
			dtemp = c*c*d[i] + s*s*d[j] - 2.0*c*s*b;
			d[j] = s*s*d[i] + c*c*d[j] +  2.0*c*s*b;
			d[i] = dtemp;
		}
	}
	
//Exit_now:

//	nrot = l;

	for(j=0; j<3; j++) {
		k = j; dtemp = d[k];
		for (i=j+1; i<4; i++) {
			if(d[i] < dtemp) {
				k = i; dtemp = d[k];
			}
		}

		if(k > j) {
			d[k] = d[j];
			d[j] = dtemp;
			for (i = 0; i <= 3; i++) {
				dtemp = v[i][k];
				v[i][k] = v[i][j];
				v[i][j] = dtemp;
			}
		}
	}
}
void quatfit(int n, vector<double> &w, vector<Xvec> &x1, vector<Xvec> &x2, double u[][3]){
	double q[4];
	double xxyx, xxyy, xxyz;
	double xyyx, xyyy, xyyz;
	double xzyx, xzyy, xzyz;
	double c[4][4], v[4][4];
	double d[4];
// generate the upper triangle of the quadratic form matrix
	xxyx = xxyy = xxyz = 0.;
	xyyx = xyyy = xyyz = 0.;
	xzyx = xzyy = xzyz = 0.;
	
	for(int i=0; i<n; i++){
		xxyx += x2[i][0] * x1[i][0] * w[i];
		xxyy += x2[i][0] * x1[i][1] * w[i];
		xxyz += x2[i][0] * x1[i][2] * w[i];
		xyyx += x2[i][1] * x1[i][0] * w[i];
		xyyy += x2[i][1] * x1[i][1] * w[i];
		xyyz += x2[i][1] * x1[i][2] * w[i];
		xzyx += x2[i][2] * x1[i][0] * w[i];
		xzyy += x2[i][2] * x1[i][1] * w[i];
		xzyz += x2[i][2] * x1[i][2] * w[i];
	}

	c[0][0] = xxyx + xyyy + xzyz;
	c[0][1] = xzyy - xyyz;
	c[1][1] = xxyx - xyyy - xzyz;
	c[0][2] = xxyz - xzyx;
	c[1][2] = xxyy + xyyx;
	c[2][2] = xyyy - xzyz - xxyx;
	c[0][3] = xyyx - xxyy;
	c[1][3] = xzyx + xxyz;
	c[2][3] = xyyz + xzyy;
	c[3][3] = xzyz - xxyx - xyyy;

// diagonalize c
	jacobi (c, d, v);

// extract the desired quaternion
	q[0] = v[0][3];
	q[1] = v[1][3];
	q[2] = v[2][3];
	q[3] = v[3][3];

// generate the rotation matrix
	u[0][0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
	u[1][0] = 2.0 * (q[1] * q[2] - q[0] * q[3]);
	u[2][0] = 2.0 * (q[1] * q[3] + q[0] * q[2]);

	u[0][1] = 2.0 * (q[2] * q[1] + q[0] * q[3]);
	u[1][1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
	u[2][1] = 2.0 * (q[2] * q[3] - q[0] * q[1]);

	u[0][2] = 2.0 *(q[3] * q[1] - q[0] * q[2]);
	u[1][2] = 2.0 * (q[3] * q[2] + q[0] * q[1]);
	u[2][2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
}
void translate(int na, vector<Xvec> &x0, vector<Xvec> &xn, double *xc){
	for(int m=0; m<3; m++) xc[m] = 0.;
	for(int i=0; i<na; i++){
		for(int m=0; m<3; m++) xc[m] += x0[i][m];
	}
	for(int m=0; m<3; m++) xc[m] /= na;
	for(int i=0; i<na; i++){
		for(int m=0; m<3; m++) xn[i][m] = x0[i][m] - xc[m];
	}
}
void translate2(int na, vector<double> &wfit, vector<Xvec> &x0, vector<Xvec> &xn, double *xc){
	double ds = 0.;
	for(int m=0; m<3; m++) xc[m] = 0.;
	for(int i=0; i<na; i++){
		for(int m=0; m<3; m++) xc[m] += x0[i][m] * wfit[i];
		ds += wfit[i];
	}
	for(int m=0; m<3; m++) xc[m] /= ds;
	for(int i=0; i<na; i++){
		for(int m=0; m<3; m++) xn[i][m] = x0[i][m] - xc[m];
	}
}
void rotmol(int na, double u[][3], vector<Xvec> &x0, vector<Xvec> &xn){
	double xt[3];
	for(int i=0; i<na; i++){
		for(int m=0; m<3; m++) xt[m] = dot_product(u[m], x0[i].getx());
		for(int m=0; m<3; m++) xn[i][m] = xt[m] + u[3][m];
	}
}

void rotmol_ideal(int nsa, double u[][3], vector<Xvec> &x0, vector<Xvec> &xn, vector<int> &start, vector<int> &mid, vector<int> &end, vector<int> &ssec){
	double xt[3];
	#pragma ivdep
	for(int i=0; i<nsa; i++){
			for(int m=0; m<3; m++) xt[m] = dot_product(u[m], x0[start[i]].getx());
			for(int m=0; m<3; m++) xn[start[i]][m] = xt[m] + u[3][m];
			for(int m=0; m<3; m++) xt[m] = dot_product(u[m], x0[mid[i]].getx());
			for(int m=0; m<3; m++) xn[mid[i]][m] = xt[m] + u[3][m];
			for(int m=0; m<3; m++) xt[m] = dot_product(u[m], x0[end[i]].getx());
			for(int m=0; m<3; m++) xn[end[i]][m] = xt[m] + u[3][m];
	}
}

inline double distance2(double *xap, double *xbp){
	double xd[3];
	for(int m=0; m<3; m++) xd[m] = xbp[m] - xap[m];
	return dot_product(xd, xd);
}
