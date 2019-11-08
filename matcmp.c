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
	
	tmp = (double)(*counts1) / tot1 - (double)(*counts2) / tot2;
	d = tmp < 0 ? -tmp : tmp;
	i = 6;
	while(--i) {
		tmp = (double)(*++counts1) / tot1 - (double)(*++counts2) / tot2;
		d += tmp < 0 ? -tmp : tmp;
	}
	
	return d;
}

double nl2cmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	int i;
	double d, tmp;
	
	tmp = (double)(*counts1) / tot1 - (double)(*counts2) / tot2;
	d = tmp * tmp;
	i = 6;
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
	
	d = pow((double)(*counts1) / tot1 - (double)(*counts2) / tot2, n);
	i = 6;
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
	
	tmp = (double) *counts1 / tot1 - (double) *counts2 / tot2;
	max = tmp < 0 ? -tmp : tmp;
	i = 6;
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
	i = 6;
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
	i = 6;
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
	i = 6;
	while(--i) {
		d += pow(abs(*++counts1 - *++counts2), n);
	}
	d = pow(d, 1.0 / n);
	
	return d < 0 ? 0 : d;
}

double linfcmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	tot1 = abs(*counts1 - *counts2);
	tot2 = 6;
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
	
	t1 = (double) *counts1 / tot1;
	t2 = (double) *counts2 / tot2;
	d = t1 < t2 ? t1 : t2;
	i = 6;
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
	i = 6;
	while(--i) {
		d += *++counts1 < *++counts2 ? *counts1 : *counts2;
	}
	d /= (tot1 + tot2);
	d =  1 - 2 * d;
	
	return d < 0 ? 0 : d;
}

double zcmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	static double alpha = 0.05;
	int i, x1, x2, max1, max2;
	
	if(counts1 == 0) {
		alpha = *((double *)(counts2));
		return 0;
	}
	
	x1 = 0;
	x2 = 0;
	max1 = *counts1;
	max2 = *counts2;
	i = 6;
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
	
	if(tot1 < (max1 << 1) && tot2 < (max2 << 1) && 
		p_chisqr(pow(tot1 - (max1 << 1), 2) / tot1) <= alpha && 
		p_chisqr(pow(tot2 - (max2 << 1), 2) / tot2) <= alpha) {
		
		return x1 == x2 ? 0 : 1;
	}
	
	return -1;
}

double chi2cmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	double d;
	
	d = (tot1 = *counts1 - *counts2) ? (tot1 * tot1 / ((double) (*counts1 + *counts2))) : 0;
	tot2 = 6;
	while(--tot2) {
		if((tot1 = *++counts1 - *++counts2)) {
			d += tot1 * tot1 / ((double) (*counts1 + *counts2));
		}
	}
	
	return sqrt(d);
}

double nchi2cmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2) {
	
	int n;
	double d, diff, t1, t2;
	
	t1 = (double) *counts1 / tot1;
	t2 = (double) *counts2 / tot2;
	d = (diff = t1 - t2) != 0 ? (diff * diff / (t1 + t2)) : 0;
	n = 6;
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
	
	int i;
	double d;
	
	tot1 = *counts1 * *counts1;
	tot2 = *counts2 * *counts2;
	d = *counts1 * *counts2;
	i = 6;
	while(--i) {
		d += *++counts1 * *++counts2;
		tot1 += *counts1 * *counts1;
		tot2 += *counts2 * *counts2;
	}
	d = 1 - d / (sqrt(tot1) * sqrt(tot2));
	
	/* take care of negative approx. error */
	return d < 0 ? 0 : d;
}

double cmpMats(MatrixCounts *mat1, NucCount *mat2, FileBuff *infile, unsigned norm, unsigned minDepth, unsigned minLength, double minCov, double (*veccmp)(short unsigned*, short unsigned*, int, int)) {
	
	/* requires mat1 to be stripped from insertions */
	
	unsigned rowNum, rowsInc, nNucs;
	short unsigned *counts1;
	double dist, d, tot1;
	
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
				if(minDepth <= (tot1 = counts1[6]) && 0 <= (d = veccmp(counts1, mat2->counts, tot1, mat2->total))) {
					dist += d;
					++rowsInc;
				}
			}
			counts1 += 7;
		}
	}
	
	if(nNucs < minLength || nNucs < minCov * rowNum) {
		return -2.0;
	} else if(rowsInc < minLength) {
		mat2->total = 0;
		return -1.0;
	}
	
	mat2->total = nNucs;
	return dist / rowsInc * norm;
}
