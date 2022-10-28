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

#define _XOPEN_SOURCE 600
#include <stdlib.h>
#include <stdio.h>
#include <sys/mman.h>
#include <sys/param.h>
#include "dat.h"
#include "pherror.h"
#include "tmp.h"

Dat * (*Dat_init)(int, int) = &DatInit;

Dat * DatInit(int M, int N) {
	
	static long unsigned type = sizeof(double);
	long unsigned Size;
	double **ptr, *src;
	float **fptr, *fsrc;
	short unsigned **sptr, *ssrc;
	unsigned char **bptr, *bsrc;
	Dat *dest;
	
	if(M < 0) {
		type = -M;
		return 0;
	}
	
	dest = smalloc(sizeof(Dat));
	dest->m = 0;
	dest->n = 0;
	dest->M= M;
	dest->N= N;
	Size = M;
	Size *= N;
	Size *= type;
	dest->mat = 0;
	dest->fmat = 0;
	dest->smat = 0;
	dest->bmat = 0;
	if(type == sizeof(double)) {
		dest->mat = smalloc(M * sizeof(double *));
		*(dest->mat) = smalloc(Size);
	} else if(type == sizeof(float)) {
		dest->fmat = smalloc(M * sizeof(float *));
		*(dest->fmat) = smalloc(Size);
	} else if(type == sizeof(short unsigned)) {
		dest->smat = smalloc(M * sizeof(short unsigned *));
		*(dest->smat) = smalloc(Size);
	} else {
		dest->bmat = smalloc(M * sizeof(unsigned char *));
		*(dest->bmat) = smalloc(Size);
	}
	dest->file = 0;
	
	/* set matrix rows */
	if(dest->mat) {
		ptr = dest->mat;
		src = *ptr;
		*ptr = src;
		while(--M) {
			*++ptr = (src += N);
		}
	} else if(dest->fmat) {
		fptr = dest->fmat;
		fsrc = *fptr;
		*fptr = fsrc;
		while(--M) {
			*++fptr = (fsrc += N);
		}
	} else if(dest->smat) {
		sptr = dest->smat;
		ssrc = *sptr;
		*sptr = ssrc;
		while(--M) {
			*++sptr = (ssrc += N);
		}
	} else {
		bptr = dest->bmat;
		bsrc = *bptr;
		*bptr = bsrc;
		while(--M) {
			*++bptr = (bsrc += N);
		}
	}
	
	return dest;
}

Dat * DatMinit(int M, int N) {
	
	static long unsigned type = sizeof(double);
	long unsigned Size;
	double **ptr, *src;
	float **fptr, *fsrc;
	short unsigned **sptr, *ssrc;
	unsigned char **bptr, *bsrc;
	FILE *tmp;
	Dat *dest;
	
	if(M < 0) {
		type = -M;
		return 0;
	}
	
	dest = smalloc(sizeof(Dat));
	dest->m = 0;
	dest->n = 0;
	dest->M= M;
	dest->N= N;
	dest->mat = 0;
	dest->fmat = 0;
	dest->smat = 0;
	dest->bmat = 0;
	if(type == sizeof(double)) {
		dest->mat = smalloc(M * sizeof(double *));
	} else if(type == sizeof(float)) {
		dest->fmat = smalloc(M * sizeof(float *));
	} else if(type == sizeof(short unsigned)) {
		dest->smat = smalloc(M * sizeof(short unsigned *));
	} else {
		dest->bmat = smalloc(M * sizeof(unsigned char *));
	}
	tmp = tmpF(0);
	dest->file = tmp;
	Size = M;
	Size *= N;
	Size *= type;
	if(fseek(tmp, Size - 1, SEEK_SET) || putc(0, tmp) == EOF) {
		ERROR();
	}
	fflush(tmp);
	sfseek(tmp, 0, SEEK_SET);
	if(dest->mat) {
		*(dest->mat) = mmap(0, Size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(tmp), 0);
		if(*(dest->mat) == MAP_FAILED) {
			fprintf(stderr, "MMAP failed:\n");
			ERROR();
		}
		//posix_madvise(*(dest->mat), Size, POSIX_MADV_SEQUENTIAL);
		
		/* set matrix rows */
		ptr = dest->mat;
		src = *ptr;
		*ptr = src;
		while(--M) {
			*++ptr = (src += N);
		}
	} else if(dest->fmat) {
		*(dest->fmat) = mmap(0, Size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(tmp), 0);
		if(*(dest->fmat) == MAP_FAILED) {
			fprintf(stderr, "MMAP failed:\n");
			ERROR();
		}
		//posix_madvise(*(dest->fmat), Size, POSIX_MADV_SEQUENTIAL);
		
		/* set matrix rows */
		fptr = dest->fmat;
		fsrc = *fptr;
		*fptr = fsrc;
		while(--M) {
			*++fptr = (fsrc += N);
		}
	} else if(dest->smat) {
		*(dest->smat) = mmap(0, Size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(tmp), 0);
		if(*(dest->smat) == MAP_FAILED) {
			fprintf(stderr, "MMAP failed:\n");
			ERROR();
		}
		//posix_madvise(*(dest->smat), Size, POSIX_MADV_SEQUENTIAL);
		
		/* set matrix rows */
		sptr = dest->smat;
		ssrc = *sptr;
		*sptr = ssrc;
		while(--M) {
			*++sptr = (ssrc += N);
		}
	} else {
		*(dest->bmat) = mmap(0, Size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(tmp), 0);
		if(*(dest->bmat) == MAP_FAILED) {
			fprintf(stderr, "MMAP failed:\n");
			ERROR();
		}
		//posix_madvise(*(dest->bmat), Size, POSIX_MADV_SEQUENTIAL);
		
		/* set matrix rows */
		bptr = dest->bmat;
		bsrc = *bptr;
		*bptr = bsrc;
		while(--M) {
			*++bptr = (bsrc += N);
		}
	}
	
	return dest;
}

void Dat_mrealloc(Dat *src, int M) {
	
	long unsigned type, Size;
	double **ptr, *mat;
	float **fptr, *fmat;
	short unsigned **sptr, *smat;
	unsigned char **bptr, *bmat;
	FILE *tmp;
	
	/* init */
	type = src->mat ? sizeof(double) : src->fmat ? sizeof(float) : src->smat ? sizeof(short unsigned) : sizeof(unsigned char);
	
	/* unmap current mapping */
	Size = src->M;
	Size *= src->N;
	Size *= type;
	if(src->mat) {
		msync(*(src->mat), Size, MS_SYNC);
		munmap(*(src->mat), Size);
	} else if(src->fmat) {
		msync(*(src->fmat), Size, MS_SYNC);
		munmap(*(src->fmat), Size);	
	} else if(src->smat) {
		msync(*(src->smat), Size, MS_SYNC);
		munmap(*(src->smat), Size);	
	} else {
		msync(*(src->bmat), Size, MS_SYNC);
		munmap(*(src->bmat), Size);	
	}
	
	/* reallocate file */
	tmp = src->file;
	Size = M;
	Size *= src->N;
	Size *= type;
	if(fseek(tmp, Size - 1, SEEK_SET) || putc(0, tmp) == EOF) {
		ERROR();
	}
	fflush(tmp);
	sfseek(tmp, 0, SEEK_SET);
	
	/* map new size */
	if(type == sizeof(double)) {
		mat = mmap(0, Size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(tmp), 0);
		if(mat == MAP_FAILED) {
			fprintf(stderr, "MMAP failed:\n");
			ERROR();
		}
		//posix_madvise(mat, Size, POSIX_MADV_SEQUENTIAL);
		ptr = realloc(src->mat, M * sizeof(double *));
		if(!ptr) {
			ERROR();
		}
		src->mat = ptr;
		*(src->mat) = mat;
		src->M = M;
		
		/* set matrix rows */
		*ptr = mat;
		while(--M) {
			*++ptr = (mat += src->N);
		}
	} else if(type == sizeof(float)) {
		fmat = mmap(0, Size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(tmp), 0);
		if(fmat == MAP_FAILED) {
			fprintf(stderr, "MMAP failed:\n");
			ERROR();
		}
		//posix_madvise(fmat, Size, POSIX_MADV_SEQUENTIAL);
		fptr = realloc(src->fmat, M * sizeof(float *));
		if(!fptr) {
			ERROR();
		}
		src->fmat = fptr;
		*(src->fmat) = fmat;
		src->M = M;
		
		/* set matrix rows */
		*fptr = fmat;
		while(--M) {
			*++fptr = (fmat += src->N);
		}
	} else if(type == sizeof(short unsigned)) {
		smat = mmap(0, Size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(tmp), 0);
		if(smat == MAP_FAILED) {
			fprintf(stderr, "MMAP failed:\n");
			ERROR();
		}
		//posix_madvise(smat, Size, POSIX_MADV_SEQUENTIAL);
		sptr = realloc(src->smat, M * sizeof(short unsigned *));
		if(!sptr) {
			ERROR();
		}
		src->smat = sptr;
		*(src->smat) = smat;
		src->M = M;
		
		/* set matrix rows */
		*sptr = smat;
		while(--M) {
			*++sptr = (smat += src->N);
		}
	} else {
		bmat = mmap(0, Size, PROT_READ | PROT_WRITE, MAP_SHARED, fileno(tmp), 0);
		if(bmat == MAP_FAILED) {
			fprintf(stderr, "MMAP failed:\n");
			ERROR();
		}
		//posix_madvise(bmat, Size, POSIX_MADV_SEQUENTIAL);
		bptr = realloc(src->bmat, M * sizeof(unsigned char *));
		if(!bptr) {
			ERROR();
		}
		src->bmat = bptr;
		*(src->bmat) = bmat;
		src->M = M;
		
		/* set matrix rows */
		*bptr = bmat;
		while(--M) {
			*++bptr = (bmat += src->N);
		}
	}
}

void Dat_realloc(Dat *src, int M) {
	
	long unsigned Size;
	double **ptr, *mat;
	float **fptr, *fmat;
	short unsigned **sptr, *smat;
	unsigned char **bptr, *bmat;
	
	if(src->file) {
		Dat_mrealloc(src, M);
		return;
	}
	
	Size = M;
	Size *= src->N;
	if(src->mat) {
		Size *= sizeof(double);
		mat = realloc(*(src->mat), Size);
		ptr = realloc(src->mat, M * sizeof(double *));
		if(!ptr || !mat) {
			ERROR();
		}
		src->mat = ptr;
		*(src->mat) = mat;
		src->M = M;
		
		/* set matrix rows */
		*ptr = mat;
		while(--M) {
			*++ptr = (mat += src->N);
		}
	} else if(src->fmat) {
		Size *= sizeof(float);
		fmat = realloc(*(src->fmat), Size);
		fptr = realloc(src->fmat, M * sizeof(float *));
		if(!fptr || !fmat) {
			ERROR();
		}
		src->fmat = fptr;
		*(src->fmat) = fmat;
		src->M = M;
		
		/* set matrix rows */
		*fptr = fmat;
		while(--M) {
			*++fptr = (fmat += src->N);
		}
	} else if(src->smat) {
		Size *= sizeof(short unsigned);
		smat = realloc(*(src->smat), Size);
		sptr = realloc(src->smat, M * sizeof(short unsigned *));
		if(!sptr || !smat) {
			ERROR();
		}
		src->smat = sptr;
		*(src->smat) = smat;
		src->M = M;
		
		/* set matrix rows */
		*sptr = smat;
		while(--M) {
			*++sptr = (smat += src->N);
		}
	} else {
		Size *= sizeof(unsigned char);
		bmat = realloc(*(src->bmat), Size);
		bptr = realloc(src->bmat, M * sizeof(unsigned char *));
		if(!bptr || !bmat) {
			ERROR();
		}
		src->bmat = bptr;
		*(src->bmat) = bmat;
		src->M = M;
		
		/* set matrix rows */
		*bptr = bmat;
		while(--M) {
			*++bptr = (bmat += src->N);
		}
	}
}

void Dat_mdestroy(Dat *src) {
	
	long unsigned size;
	
	if(src) {
		size = src->M;
		size *= src->N;
		if(src->mat) {
			size *= sizeof(double);
			msync(*(src->mat), size, MS_SYNC);
			munmap(*(src->mat), size);
			free(src->mat);
		} else if(src->fmat) {
			size *= sizeof(float);
			msync(*(src->fmat), size, MS_SYNC);
			munmap(*(src->fmat), size);
			free(src->fmat);
		} else if(src->smat) {
			size *= sizeof(short unsigned);
			msync(*(src->smat), size, MS_SYNC);
			munmap(*(src->smat), size);
			free(src->smat);
		} else if(src->bmat) {
			size *= sizeof(unsigned char);
			msync(*(src->bmat), size, MS_SYNC);
			munmap(*(src->bmat), size);
			free(src->bmat);
		}
		fclose(src->file);
		free(src);
	}
}

void Dat_destroy(Dat *src) {
	
	if(src) {
		if(src->file) {
			Dat_mdestroy(src);
			return;
		} else if(src->mat) {
			free(*(src->mat));
			free(src->mat);
		} else if(src->fmat) {
			free(*(src->fmat));
			free(src->fmat);
		} else if(src->smat) {
			free(*(src->smat));
			free(src->smat);
		} else if(src->bmat) {
			free(*(src->bmat));
			free(src->bmat);
		}
		free(src);
	}
}
