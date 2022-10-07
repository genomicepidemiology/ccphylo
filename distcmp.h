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

extern double (*distcmp_d)(double *v1, double *v2, long n);
extern double (*distcmp_f)(float *v1, float *v2, long n);
extern double (*distcmp_s)(short unsigned *v1, short unsigned *v2, long n);
extern double (*distcmp_b)(unsigned char *v1, unsigned char *v2, long n);

double l1cmp_d(double *v1, double *v2, long n);
double l1cmp_f(float *v1, float *v2, long n);
double l1cmp_s(short unsigned *v1, short unsigned *v2, long n);
double l1cmp_b(unsigned char *v1, unsigned char *v2, long n);
double l2cmp_d(double *v1, double *v2, long n);
double l2cmp_f(float *v1, float *v2, long n);
double l2cmp_s(short unsigned *v1, short unsigned *v2, long n);
double l2cmp_b(unsigned char *v1, unsigned char *v2, long n);
double lncmp_d(double *v1, double *v2, long n);
double lncmp_f(float *v1, float *v2, long n);
double lncmp_s(short unsigned *v1, short unsigned *v2, long n);
double lncmp_b(unsigned char *v1, unsigned char *v2, long n);
double linfcmp_d(double *v1, double *v2, long n);
double linfcmp_f(float *v1, float *v2, long n);
double linfcmp_s(short unsigned *v1, short unsigned *v2, long n);
double linfcmp_b(unsigned char *v1, unsigned char *v2, long n);
double bccmp_d(double *v1, double *v2, long n);
double bccmp_f(float *v1, float *v2, long n);
double bccmp_s(short unsigned *v1, short unsigned *v2, long n);
double bccmp_b(unsigned char *v1, unsigned char *v2, long n);
double chi2cmp_d(double *v1, double *v2, long n);
double chi2cmp_f(float *v1, float *v2, long n);
double chi2cmp_s(short unsigned *v1, short unsigned *v2, long n);
double chi2cmp_b(unsigned char *v1, unsigned char *v2, long n);
double coscmp_d(double *v1, double *v2, long n);
double coscmp_f(float *v1, float *v2, long n);
double coscmp_s(short unsigned *v1, short unsigned *v2, long n);
double coscmp_b(unsigned char *v1, unsigned char *v2, long n);
double pearcmp_d(double *v1, double *v2, long n);
double pearcmp_f(float *v1, float *v2, long n);
double pearcmp_s(short unsigned *v1, short unsigned *v2, long n);
double pearcmp_b(unsigned char *v1, unsigned char *v2, long n);
