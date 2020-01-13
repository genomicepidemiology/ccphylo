/* Philip T.L.C. Clausen Jan 2017 plan@dtu.dk */

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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "filebuff.h"
#include "matparse.h"
#include "stdstat.h"

void stripMat(MatrixCounts *mat1) {
	
	unsigned char *ref, *validRef;
	unsigned i, j, len;
	short unsigned *valid, *ptr;
	
	ref = mat1->refs;
	i = mat1->len + 1;
	len = 0;
	while(i && *ref != '-') {
		--i;
		++ref;
		++len;
	}
	
	if(i) {
		validRef = ref;
		valid = mat1->counts + (7 * len);
		ptr = --valid;
		while(--i) {
			if(*ref != '-') {
				*validRef++ = *ref++;
				j = 8;
				while(--j) {
					*++valid = *++ptr;
				}
				++len;
			} else {
				++ref;
				ptr += 7;
			}
		}
	}
	mat1->len = len;
}

double nl1cmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	int i;
	double d, tmp;
	
	tot1 -= counts1[5];
	tot2 -= counts2[5];
	tmp = (double)(*counts1) / tot1 - (double)(*counts2) / tot2;
	d = tmp < 0 ? -tmp : tmp;
	i = 5;
	while(--i) {
		tmp = (double)(*++counts1) / tot1 - (double)(*++counts2) / tot2;
		d += tmp < 0 ? -tmp : tmp;
	}
	
	return d;
}

double nl2cmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	int i;
	double d, tmp;
	
	tot1 -= counts1[5];
	tot2 -= counts2[5];
	tmp = (double)(*counts1) / tot1 - (double)(*counts2) / tot2;
	d = tmp * tmp;
	i = 5;
	while(--i) {
		tmp = (double)(*++counts1) / tot1 - (double)(*++counts2) / tot2;
		d += tmp * tmp;
	}
	
	return sqrt(d);
}

double nlncmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	static unsigned n = 0;
	int i;
	double d, tmp;
	
	if(counts1 == 0) {
		n = *((unsigned *)(counts2));
		return 0;
	}
	
	tot1 -= counts1[5];
	tot2 -= counts2[5];
	d = pow((double)(*counts1) / tot1 - (double)(*counts2) / tot2, n);
	i = 5;
	while(--i) {
		tmp = (double)(*++counts1) / tot1 - (double)(*++counts2) / tot2;
		tmp = tmp < 0 ? -tmp : tmp;
		d += pow(tmp, n);
	}
	d = pow(d, 1.0 / n);
	
	return d < 0 ? 0 : d;
}

double nlinfcmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	int i;
	double max, tmp;
	
	tot1 -= counts1[5];
	tot2 -= counts2[5];
	tmp = (double) *counts1 / tot1 - (double) *counts2 / tot2;
	max = tmp < 0 ? -tmp : tmp;
	i = 5;
	while(--i) {
		tmp = (double) *counts1 / tot1 - (double) *counts2 / tot2;
		tmp = tmp < 0 ? -tmp : tmp;
		if(max < tmp) {
			max = tmp;
		}
	}
	
	return max;
}

double l1cmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	int i;
	
	tot1 = *counts1 - *counts2;
	tot2 = tot1 < 0 ? -tot1 : tot1;
	i = 5;
	while(--i) {
		tot1 = *++counts1 - *++counts2;
		tot2 += tot1 < 0 ? -tot1 : tot1;
	}
	
	return tot2;
}

double l2cmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	int i;
	
	tot1 = *counts1 - *counts2;
	tot2 = tot1 * tot1;
	i = 5;
	while(--i) {
		tot1 = *++counts1 - *++counts2;
		tot2 += tot1 * tot1;
	}
	
	return sqrt(tot2);
}

double lncmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	static unsigned n = 0;
	int i;
	double d;
	
	if(counts1 == 0) {
		n = *((unsigned *)(counts2));
		return 0;
	}
	
	d = pow(abs(*counts1 - *counts2), n);
	i = 5;
	while(--i) {
		d += pow(abs(*++counts1 - *++counts2), n);
	}
	d = pow(d, 1.0 / n);
	
	return d < 0 ? 0 : d;
}

double linfcmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	tot1 = abs(*counts1 - *counts2);
	tot2 = 5;
	while(--tot2) {
		if(tot1 < abs(*++counts1 - *++counts2)) {
			tot1 = abs(*counts1 - *counts2);
		}
	}
	
	return tot1;
}

double nbccmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	int i;
	double d, t1, t2;
	
	tot1 -= counts1[5];
	tot2 -= counts2[5];
	t1 = (double) *counts1 / tot1;
	t2 = (double) *counts2 / tot2;
	d = t1 < t2 ? t1 : t2;
	i = 5;
	while(--i) {
		t1 = (double) *++counts1 / tot1;
		t2 = (double) *++counts2 / tot2;
		d += t1 < t2 ? t1 : t2;
	}
	d = 1 - d;
	
	return d < 0 ? 0 : d;
}

double bccmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	int i;
	double d;
	
	d = *counts1 < *counts2 ? *counts1 : *counts2;
	i = 5;
	while(--i) {
		d += *++counts1 < *++counts2 ? *counts1 : *counts2;
	}
	d /= (tot1 - *++counts1 + tot2 - *++counts2);
	d =  1 - 2 * d;
	
	return d < 0 ? 0 : d;
}

double nccmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	int i;
	double d, t1, t2, T;
	
	tot1 -= counts1[5];
	tot2 -= counts2[5];
	t1 = (double) *counts1 / tot1;
	t2 = (double) *counts2 / tot2;
	d = 1;
	if(t1 < t2) {
		d = t1;
		T = t2;
	} else {
		d = t2;
		T = t1;
	}
	i = 5;
	while(--i) {
		t1 = (double) *++counts1 / tot1;
		t2 = (double) *++counts2 / tot2;
		T = 1;
		if(t1 < t2) {
			d += t1;
			T += t2;
		} else {
			d += t2;
			T += t1;
		}
	}
	d = 1 - d / T;
	
	return d < 0 ? 0 : d;
}

double ccmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	int i;
	double d;
	
	if(*counts1 < *counts2) {
		d = *counts1;
		tot1 = *counts2;
	} else {
		d = *counts2;
		tot1 = *counts1;
	}
	i = 5;
	while(--i) {
		if(*++counts1 < *++counts2) {
			d += *counts1;
			tot1 += *counts2;
		} else {
			d += *counts2;
			tot1 += *counts1;
		}
	}
	if(!tot1) {
		return -1;
	}
	d = 1 - d / tot1;
	
	return d < 0 ? 0 : d;
}

double zcmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	static double alpha = 0.05;
	int i, x1, x2, max1, max2;
	
	if(counts1 == 0) {
		alpha = *((double *)(counts2));
		return 0;
	}
	
	i = 5;
	x1 = 0;
	x2 = 0;
	max1 = *counts1;
	max2 = *counts2;
	while(--i) {
		if(max1 < *++counts1) {
			max1 = *counts1;
			x1 = i;
		}
		if(max2 < *++counts2) {
			max2 = *counts2;
			x2 = i;
		}
	}
	
	x1 = p_chisqr(pow(tot1 - (max1 << 1), 2) / tot1) <= alpha && tot1 < (max1 << 1);
	x2 = p_chisqr(pow(tot2 - (max2 << 1), 2) / tot2) <= alpha && tot1 < (max1 << 1);
	if(x1 && x2) {
		return x1 == x2 ? 0 : 1;
	}
	
	return -1;
}

double pcmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	double d, T;
	
	d = (T = *counts1 - *counts2) ? T * T / (*counts1 + *counts2) : 0;
	tot1 = 5;
	while(--tot1) {
		if((T = *++counts1 - *++counts2)) {
			d += T * T / (*counts1 + *counts2);
		}
	}
	
	return 1 - p_chisqr(d);
}

double npcmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	int n;
	double d, diff, t1, t2;
	
	tot1 -= counts1[5];
	tot2 -= counts2[5];
	t1 = (double) *counts1 / tot1;
	t2 = (double) *counts2 / tot2;
	d = (diff = t1 - t2) ? diff * diff / (t1 + t2) : 0;
	n = 5;
	while(--n) {
		t1 = (double) *++counts1 / tot1;
		t2 = (double) *++counts2 / tot2;
		if((diff = t1 - t2)) {
			d += diff * diff / (t1 + t2);
		}
	}
	
	return 1 - p_chisqr(d);
}

double chi2cmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	double d, T;
	
	d = (T = *counts1 - *counts2) ? (T * T / (*counts1 + *counts2)) : 0;
	tot2 = 5;
	while(--tot2) {
		if((T = *++counts1 - *++counts2)) {
			d += T * T / (*counts1 + *counts2);
		}
	}
	
	return sqrt(d);
}

double nchi2cmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	int n;
	double d, diff, t1, t2;
	
	tot1 -= counts1[5];
	tot2 -= counts2[5];
	t1 = (double) *counts1 / tot1;
	t2 = (double) *counts2 / tot2;
	d = (diff = t1 - t2) != 0 ? (diff * diff / (t1 + t2)) : 0;
	n = 5;
	while(--n) {
		t1 = (double) *++counts1 / tot1;
		t2 = (double) *++counts2 / tot2;
		if((diff = t1 - t2)) {
			d += diff * diff / (t1 + t2);
		}
	}
	
	return sqrt(d);
}

double coscmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	unsigned i;
	long unsigned c1, c2;
	double d;
	
	tot1 = *counts1;
	tot2 = *counts2;
	c1 = tot1 * tot1;
	c2 = tot2 * tot2;
	d = tot1 * tot2;
	i = 5;
	while(--i) {
		tot1 = *++counts1;
		tot2 = *++counts2;
		d += tot1 * tot2;
		c1 += tot1 * tot1;
		c2 += tot2 * tot2;
	}
	if(!c1 || !c2) {
		return -1;
	}
	d = 1 - d / (sqrt(c1) * sqrt(c2));
	
	/* take care of negative approx. error */
	return d < 0 ? 0 : d;
}

double cmpMats(MatrixCounts *mat1, NucCount *mat2, FileBuff *infile, unsigned norm, unsigned minDepth, unsigned minLength, double minCov, double (*veccmp)(short unsigned*, short unsigned*, int, int)) {
	
	/* requires mat1 to be stripped from insertions */
	
	unsigned rowNum, rowsInc, nNucs, tot1;
	short unsigned *counts1;
	double dist, d;
	
	/* check if templates matches */
	if(strcmp((char *) mat1->name->seq, (char *) mat2->name) != 0) {
		return -2.0;
	}
	
	dist = 0;
	rowNum = 0;
	rowsInc = 0;
	nNucs = 0;
	counts1 = mat1->counts;
	while(FileBuffGetRow(infile, mat2) && mat2->ref) {
		if(mat2->ref != '-') {
			/* validate position */
			if(mat1->len < ++rowNum) {
				return -1;
			}
			if(minDepth <= mat2->total) {
				++nNucs;
				tot1 = *((unsigned *)(counts1 + 6));
				if(minDepth <= tot1 && 0 <= (d = veccmp(counts1, mat2->counts, tot1, mat2->total))) {
					dist += d;
					++rowsInc;
				}
			}
			counts1 += 8;
		}
	}
	
	if(nNucs < minLength || nNucs < minCov * rowNum) {
		return -2.0;
	} else if(rowsInc < minLength || rowsInc < minCov * rowNum) {
		mat2->total = 0;
		return -1.0;
	}
	
	mat2->total = rowsInc;
	
	return norm ? dist / rowsInc * norm : dist;
}
