/* Philip T.L.C. Clausen Jul 2021 plan@dtu.dk */

/*
 * Copyright (c) 2017, Philip Clausen, Technical University of Denmark
 * All rights reserved.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *		http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include <errno.h>
#include <math.h>
#include "bytescale.h"
#include "distcmp.h"

double (*distcmp_d)(double *v1, double *v2, long n) = &coscmp_d;
double (*distcmp_f)(float *v1, float *v2, long n) = &coscmp_f;
double (*distcmp_s)(short unsigned *v1, short unsigned *v2, long n) = &coscmp_s;
double (*distcmp_b)(unsigned char *v1, unsigned char *v2, long n) = &coscmp_b;

double l1cmp_d(double *v1, double *v2, long n) {
	
	double d, tmp;
	
	tmp = *v1 - *v2;
	d = tmp < 0 ? -tmp : tmp;
	while(--n) {
		tmp = *++v1 - *++v2;
		d += tmp < 0 ? -tmp : tmp;
	}
	
	return d;
}

double l1cmp_f(float *v1, float *v2, long n) {
	
	double d, tmp;
	
	tmp = *v1 - *v2;
	d = tmp < 0 ? -tmp : tmp;
	while(--n) {
		tmp = *++v1 - *++v2;
		d += tmp < 0 ? -tmp : tmp;
	}
	
	return d;
}

double l1cmp_s(short unsigned *v1, short unsigned *v2, long n) {
	
	double d, tmp;
	
	tmp = *v1 - *v2;
	d = tmp < 0 ? -tmp : tmp;
	while(--n) {
		tmp = *++v1 - *++v2;
		d += tmp < 0 ? -tmp : tmp;
	}
	
	return uctod(d);
}

double l1cmp_b(unsigned char *v1, unsigned char *v2, long n) {
	
	double d, tmp;
	
	tmp = *v1 - *v2;
	d = tmp < 0 ? -tmp : tmp;
	while(--n) {
		tmp = *++v1 - *++v2;
		d += tmp < 0 ? -tmp : tmp;
	}
	
	return uctod(d);
}

double l2cmp_d(double *v1, double *v2, long n) {
	
	double d, tmp;
	
	tmp = *v1 - *v2;
	d = tmp * tmp;;
	while(--n) {
		tmp = *++v1 - *++v2;
		d += tmp * tmp;
	}
	
	return sqrt(d);
}

double l2cmp_f(float *v1, float *v2, long n) {
	
	double d, tmp;
	
	tmp = *v1 - *v2;
	d = tmp * tmp;;
	while(--n) {
		tmp = *++v1 - *++v2;
		d += tmp * tmp;
	}
	
	return sqrt(d);
}

double l2cmp_s(short unsigned *v1, short unsigned *v2, long n) {
	
	double d, tmp;
	
	tmp = uctod(*v1 - *v2);
	d = tmp * tmp;;
	while(--n) {
		tmp = uctod(*++v1 - *++v2);
		d += tmp * tmp;
	}
	
	return sqrt(d);
}

double l2cmp_b(unsigned char *v1, unsigned char *v2, long n) {
	
	double d, tmp;
	
	tmp = uctod(*v1 - *v2);
	d = tmp * tmp;;
	while(--n) {
		tmp = uctod(*++v1 - *++v2);
		d += tmp * tmp;
	}
	
	return sqrt(d);
}

double lncmp_d(double *v1, double *v2, long n) {
	
	static double exponent = 0;
	double d, tmp;
	
	if(v1 == 0) {
		exponent = *v2;
		return 0;
	}
	
	tmp = *v1 - *v2;
	d = pow(tmp < 0 ? -tmp : tmp, exponent);
	while(--n) {
		tmp = *++v1 - *++v2;
		d += pow(tmp < 0 ? -tmp : tmp, exponent);
	}
	d = pow(d, 1.0 / exponent);
	
	return d < 0 ? 0 : d;
}

double lncmp_f(float *v1, float *v2, long n) {
	
	static double exponent = 0;
	double d, tmp;
	
	if(v1 == 0) {
		exponent = *((double *)(v2));
		return 0;
	}
	
	tmp = *v1 - *v2;
	d = pow(tmp < 0 ? -tmp : tmp, exponent);
	while(--n) {
		tmp = *++v1 - *++v2;
		d += pow(tmp < 0 ? -tmp : tmp, exponent);
	}
	d = pow(d, 1.0 / exponent);
	
	return d < 0 ? 0 : d;
}

double lncmp_b(unsigned char *v1, unsigned char *v2, long n) {
	
	static double exponent = 0;
	double d, tmp;
	
	if(v1 == 0) {
		exponent = *((double *)(v2));
		return 0;
	}
	
	tmp = uctod(*v1 - *v2);
	d = pow(tmp < 0 ? -tmp : tmp, exponent);
	while(--n) {
		tmp = uctod(*++v1 - *++v2);
		d += pow(tmp < 0 ? -tmp : tmp, exponent);
	}
	d = pow(d, 1.0 / exponent);
	
	return d < 0 ? 0 : d;
}

double lncmp_s(short unsigned *v1, short unsigned *v2, long n) {
	
	static double exponent = 0;
	double d, tmp;
	
	if(v1 == 0) {
		exponent = *((double *)(v2));
		return 0;
	}
	
	tmp = uctod(*v1 - *v2);
	d = pow(tmp < 0 ? -tmp : tmp, exponent);
	while(--n) {
		tmp = uctod(*++v1 - *++v2);
		d += pow(tmp < 0 ? -tmp : tmp, exponent);
	}
	d = pow(d, 1.0 / exponent);
	
	return d < 0 ? 0 : d;
}

double linfcmp_d(double *v1, double *v2, long n) {
	
	double d, tmp;
	
	tmp = *v1 - *v2;
	d = tmp < 0 ? -tmp : tmp;
	while(--n) {
		tmp = *++v1 - *++v2;
		if(d < tmp) {
			d = tmp;
		} else if(tmp < 0 && d < -tmp) {
			d = -tmp;
		}
	}
	
	return d;
}

double linfcmp_f(float *v1, float *v2, long n) {
	
	double d, tmp;
	
	tmp = *v1 - *v2;
	d = tmp < 0 ? -tmp : tmp;
	while(--n) {
		tmp = *++v1 - *++v2;
		if(d < tmp) {
			d = tmp;
		} else if(tmp < 0 && d < -tmp) {
			d = -tmp;
		}
	}
	
	return d;
}

double linfcmp_s(short unsigned *v1, short unsigned *v2, long n) {
	
	unsigned char d, tmp;
	
	tmp = *v1 - *v2;
	d = tmp < 0 ? -tmp : tmp;
	while(--n) {
		tmp = *++v1 - *++v2;
		if(d < tmp) {
			d = tmp;
		} else if(tmp < 0 && d < -tmp) {
			d = -tmp;
		}
	}
	
	return uctod(d);
}

double linfcmp_b(unsigned char *v1, unsigned char *v2, long n) {
	
	unsigned char d, tmp;
	
	tmp = *v1 - *v2;
	d = tmp < 0 ? -tmp : tmp;
	while(--n) {
		tmp = *++v1 - *++v2;
		if(d < tmp) {
			d = tmp;
		} else if(tmp < 0 && d < -tmp) {
			d = -tmp;
		}
	}
	
	return uctod(d);
}

double bccmp_d(double *v1, double *v2, long n) {
	
	double d, s;
	
	d = *v1 < *v2 ? *v1 : *v2;
	s = *v1 + *v2;
	while(--n) {
		d += *++v1 < *++v2 ? *v1 : *v2;
		s += *v1 + *v2;
	}
	d = 1 - 2 * d / s;
	
	return d < 0 ? 0 : d;
}

double bccmp_f(float *v1, float *v2, long n) {
	
	double d, s;
	
	d = *v1 < *v2 ? *v1 : *v2;
	s = *v1 + *v2;
	while(--n) {
		d += *++v1 < *++v2 ? *v1 : *v2;
		s += *v1 + *v2;
	}
	d = 1 - 2 * d / s;
	
	return d < 0 ? 0 : d;
}

double bccmp_s(short unsigned *v1, short unsigned *v2, long n) {
	
	int d, s;
	
	d = *v1 < *v2 ? *v1 : *v2;
	s = (int)(*v1) + *v2;
	while(--n) {
		d += *++v1 < *++v2 ? *v1 : *v2;
		s += (int)(*v1) + *v2;
	}
	d = 1 - 2 * (double)(d) / s;
	
	return d < 0 ? 0 : d;
}

double bccmp_b(unsigned char *v1, unsigned char *v2, long n) {
	
	int d, s;
	
	d = *v1 < *v2 ? *v1 : *v2;
	s = (int)(*v1) + *v2;
	while(--n) {
		d += *++v1 < *++v2 ? *v1 : *v2;
		s += (int)(*v1) + *v2;
	}
	d = 1 - 2 * (double)(d) / s;
	
	return d < 0 ? 0 : d;
}

double chi2cmp_d(double *v1, double *v2, long n) {
	
	double d, T;
	
	d = (T = *v1 - *v2) ? (T * T / (*v1 + *v2)) : 0;
	while(--n) {
		if((T = *++v1 - *++v2)) {
			d += T * T / (*v1 + *v2);
		}
	}
	
	return sqrt(d);
}

double chi2cmp_f(float *v1, float *v2, long n) {
	
	double d, T;
	
	d = (T = *v1 - *v2) ? (T * T / (*v1 + *v2)) : 0;
	while(--n) {
		if((T = *++v1 - *++v2)) {
			d += T * T / (*v1 + *v2);
		}
	}
	
	return sqrt(d);
}

double chi2cmp_s(short unsigned *v1, short unsigned *v2, long n) {
	
	double d, T;
	
	d = (T = *v1 - *v2) ? (T * T / ((int)(*v1) + *v2)) : 0;
	while(--n) {
		if((T = *++v1 - *++v2)) {
			d += T * T / ((int)(*v1) + *v2);
		}
	}
	
	return sqrt(d);
}

double chi2cmp_b(unsigned char *v1, unsigned char *v2, long n) {
	
	double d, T;
	
	d = (T = *v1 - *v2) ? (T * T / ((int)(*v1) + *v2)) : 0;
	while(--n) {
		if((T = *++v1 - *++v2)) {
			d += T * T / ((int)(*v1) + *v2);
		}
	}
	
	return sqrt(d);
}

double coscmp_d(double *v1, double *v2, long n) {
	
	double d, c1, c2;
	
	d = *v1 * *v2;
	c1 = *v1 * *v1;
	c2 = *v2 * *v2;
	while(--n) {
		d += *++v1 * *++v2;
		c1 += *v1 * *v1;
		c2 += *v2 * *v2;
	}
	if(!c1 || !c2) {
		errno = EDOM;
		return -1;
	}
	d = 1 - d / sqrt(c1 * c2);
	
	/* take care of negative approx. error */
	return d < 0 ? 0 : d;
}

double coscmp_f(float *v1, float *v2, long n) {
	
	double d, c1, c2;
	
	d = *v1 * *v2;
	c1 = *v1 * *v1;
	c2 = *v2 * *v2;
	while(--n) {
		d += *++v1 * *++v2;
		c1 += *v1 * *v1;
		c2 += *v2 * *v2;
	}
	if(!c1 || !c2) {
		errno = EDOM;
		return -1;
	}
	d = 1 - d / sqrt(c1 * c2);
	
	/* take care of negative approx. error */
	return d < 0 ? 0 : d;
}

double coscmp_s(short unsigned *v1, short unsigned *v2, long n) {
	
	double d, c1, c2, t1, t2;
	
	t1 = uctod(*v1);
	t2 = uctod(*v2);
	d = t1 * t2;
	c1 = t1 * t1;
	c2 = t2 * t2;
	while(--n) {
		t1 = uctod(*++v1);
		t2 = uctod(*++v2);
		d += t1 * t2;
		c1 += t1 * t1;
		c2 += t2 * t2;
	}
	if(!c1 || !c2) {
		errno = EDOM;
		return -1;
	}
	d = 1 - d / sqrt(c1 * c2);
	
	/* take care of negative approx. error */
	return d < 0 ? 0 : d;
}

double coscmp_b(unsigned char *v1, unsigned char *v2, long n) {
	
	double d, c1, c2, t1, t2;
	
	t1 = uctod(*v1);
	t2 = uctod(*v2);
	d = t1 * t2;
	c1 = t1 * t1;
	c2 = t2 * t2;
	while(--n) {
		t1 = uctod(*++v1);
		t2 = uctod(*++v2);
		d += t1 * t2;
		c1 += t1 * t1;
		c2 += t2 * t2;
	}
	if(!c1 || !c2) {
		errno = EDOM;
		return -1;
	}
	d = 1 - d / sqrt(c1 * c2);
	
	/* take care of negative approx. error */
	return d < 0 ? 0 : d;
}

double pearcmp_d(double *v1, double *v2, long n) {
	
	long i;
	double e1, e2, v11, v12, v22;
	
	/* v's are converted to ranks for spearman */
	/*
	Pab = Cov(a,b) / sqrt(Var(a) * Var(b))
	= (E(ab) - E(a) E(b)) / sqrt((E(aa) - E(a) E(a)) (E(bb) - E(b) E(b)))
	= (Eab / n - EaEb / n^2) / sqrt((Eaa / n - EaEa / n^2) (Ebb / n - EbEb / n^2))
	= (Eab - EaEb / n) / sqrt((Eaa - EaEa / n) (Ebb - EbEb / n))
	*/
	
	/* get sums */
	e1 = *v1;
	e2 = *v2;
	v11 = *v1 * *v1;
	v12 = *v1 * *v2;
	v22 = *v2 * *v2;
	i = n;
	while(--i) {
		e1 += *++v1;
		e2 += *++v2;
		v11 += *v1 * *v1;
		v12 += *v1 * *v2;
		v22 += *v2 * *v2;
	}
	
	/* get variances */
	v11 -= e1 * e1 / n;
	v12 -= e1 * e2 / n;
	v22 -= e2 * e2 / n;
	
	/* check validity */
	if(!v11 || !v22) {
		errno = EDOM;
		return 0;
	}
	
	/* corr */
	return v12 / sqrt(v11 * v22);
}

double pearcmp_f(float *v1, float *v2, long n) {
	
	long i;
	double e1, e2, v11, v12, v22;
	
	/* get sums */
	e1 = *v1;
	e2 = *v2;
	v11 = *v1 * *v1;
	v12 = *v1 * *v2;
	v22 = *v2 * *v2;
	i = n;
	while(--i) {
		e1 += *++v1;
		e2 += *++v2;
		v11 += *v1 * *v1;
		v12 += *v1 * *v2;
		v22 += *v2 * *v2;
	}
	
	/* get variances */
	v11 -= e1 * e1 / n;
	v12 -= e1 * e2 / n;
	v22 -= e2 * e2 / n;
	
	/* check validity */
	if(!v11 || !v22) {
		errno = EDOM;
		return 0;
	}
	
	/* corr */
	return v12 / sqrt(v11 * v22);
}

double pearcmp_s(short unsigned *v1, short unsigned *v2, long n) {
	
	int t1, t2;
	long i;
	double e1, e2, v11, v12, v22;
	
	/* get sums */
	t1 = *v1;
	t2 = *v2;
	e1 = t1;
	e2 = t2;
	v11 = t1 * t1;
	v12 = t1 * t2;
	v22 = t2 * t2;
	i = n;
	while(--i) {
		t1 = *++v1;
		t2 = *++v2;
		e1 += t1;
		e2 += t2;
		v11 += t1 * t1;
		v12 += t1 * t2;
		v22 += t2 * t2;
	}
	
	/* convert to doubles */
	e1 = uctod(e1);
	e2 = uctod(e2);
	v11 = uctod(v11);
	v12 = uctod(v12);
	v22 = uctod(v22);
	
	/* get variances */
	v11 -= e1 * e1 / n;
	v12 -= e1 * e2 / n;
	v22 -= e2 * e2 / n;
	
	/* check validity */
	if(!v11 || !v22) {
		errno = EDOM;
		return 0;
	}
	
	/* corr */
	return v12 / sqrt(v11 * v22);
}

double pearcmp_b(unsigned char *v1, unsigned char *v2, long n) {
	
	int t1, t2;
	long i;
	double e1, e2, v11, v12, v22;
	
	/* get sums */
	t1 = *v1;
	t2 = *v2;
	e1 = t1;
	e2 = t2;
	v11 = t1 * t1;
	v12 = t1 * t2;
	v22 = t2 * t2;
	i = n;
	while(--i) {
		t1 = *++v1;
		t2 = *++v2;
		e1 += t1;
		e2 += t2;
		v11 += t1 * t1;
		v12 += t1 * t2;
		v22 += t2 * t2;
	}
	
	/* convert to doubles */
	e1 = uctod(e1);
	e2 = uctod(e2);
	v11 = uctod(v11);
	v12 = uctod(v12);
	v22 = uctod(v22);
	
	/* get variances */
	v11 -= e1 * e1 / n;
	v12 -= e1 * e2 / n;
	v22 -= e2 * e2 / n;
	
	/* check validity */
	if(!v11 || !v22) {
		errno = EDOM;
		return 0;
	}
	
	/* corr */
	return v12 / sqrt(v11 * v22);
}
