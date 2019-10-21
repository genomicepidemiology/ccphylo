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

#include <stdlib.h>
#include "matrix.h"
#include "pherror.h"

Matrix * ltdMatrix_init(unsigned size) {
	
	int i;
	double **ptr, *src;
	Matrix *dest;
	
	dest = smalloc(sizeof(Matrix));
	dest->n = 0;
	dest->size = size;
	dest->mat = smalloc(size * sizeof(double *));
	*(dest->mat) = smalloc(size * size * sizeof(double) / 2);
	
	/* set matrix rows */
	ptr = dest->mat;
	src = *ptr;
	i = 0;
	*ptr++ = src;
	while(--size) {
		*ptr++ = src + i;
		src += i++;
	}
	
	return dest;
}

void ltdMatrix_realloc(Matrix *src, unsigned newsize) {
	
	int i, size;
	double **ptr, *mat;
	
	*(src->mat) = realloc(*(src->mat), newsize * newsize * sizeof(double) / 2);
	src->mat = realloc(src->mat, newsize * sizeof(double *));
	if(!src->mat || !*(src->mat)) {
		ERROR();
	}
	
	/* set matrix rows */
	ptr = src->mat;
	mat = *ptr;
	size = src->size;
	i = size * (size - 1) / 2;
	ptr += size;
	*ptr++ = mat;
	size = newsize - size;
	while(--size) {
		*ptr++ = mat + i;
		mat += i++;
	}
	src->size = newsize;
}

void Matrix_destroy(Matrix *src) {
	
	free(src->mat);
	free(src);
}

void ltdMatrix_popArrange(Matrix *mat, unsigned pos) {
	
	int i, n;
	double *dest, *src;
	
	n = --mat->n;
	if(pos != n) {
		/* row to be emptied */
		dest = mat->mat[pos];
		/* row to be moved up */
		src = mat->mat[n];
		
		/* copy last row into "first" row */
		i = pos + 1;
		while(--i) {
			*dest++ = *src++;
		}
		
		/* tilt remaining part of last row into column "pos" */
		i = pos;
		++src;
		while(++i < n) {
			mat->mat[i][pos] = *src++;
		}
	}
}

int ltdMatrix_add(Matrix *src) {
	
	int i;
	double *ptr;
	
	/* realloc */
	if(++src->n == src->size) {
		ltdMatrix_realloc(src, src->size << 1);
	}
	
	/* init new row */
	i = src->n;
	ptr = src->mat[i - 1] - 1;
	while(--i) {
		*++ptr = 0;
	}
	
	return src->n - 1;
}
