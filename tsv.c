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

#include <stdlib.h>
#include <stdio.h>
#include "bytescale.h"
#include "dat.h"
#include "filebuff.h"
#include "jobs.h"
#include "pherror.h"
#include "stdstat.h"
#include "tsv.h"

Dat * loadTsv(FileBuff *infile, unsigned char sep) {
	
	int avail, size, N, n, (*buffFileBuff)(FileBuff *);
	char *msg, strbuff[32];
	unsigned char c, stop, *buff, *num, *bmat;
	double *mat;
	float *fmat;
	short unsigned *smat;
	Dat *dest;
	
	/* init */
	avail = infile->bytes;
	buff = infile->next;
	buffFileBuff = infile->buffFileBuff;
	if(avail == 0) {
		if((avail = buffFileBuff(infile)) == 0) {
			return 0;
		}
		buff = infile->buffer;
	}
	
	/* skip header and get dimensions */
	do {
		N = 1;
		while(*buff++ != '\n') {
			if(--avail == 0) {
				if((avail = buffFileBuff(infile)) == 0) {
					return 0;
				}
				buff = infile->buffer;
			}
			if(*buff == sep) {
				++N;
			}
		}
		if(--avail == 0) {
			if((avail = buffFileBuff(infile)) == 0) {
				return 0;
			}
			buff = infile->buffer;
		}
	} while(*buff == '#');
	
	/* allocate matrix */
	dest = Dat_init(1048576, N);
	dest->n = N;
	
	/* load rows */
	mat = 0;
	fmat = 0;
	smat = 0;
	bmat = 0;
	if(dest->mat) {
		mat = *(dest->mat);
	} else if(dest->fmat) {
		fmat = *(dest->fmat);
	} else if(dest->smat) {
		smat = *(dest->smat);
	} else {
		bmat = *(dest->bmat);
	}
	while(avail) {
		/* load row */
		n = N;
		while(n--) {
			stop = n ? sep : '\n';
			num = (unsigned char *) strbuff;
			size = 32;
			while((c = *buff++) != stop) {
				*num++ = c;
				if(--size == 0) {
					fprintf(stderr, "Malformatted entry at pos:\t(%d,%d) %s\n", dest->m, N - n, strbuff);
					exit(1);
				}
				if(--avail == 0) {
					if((avail = buffFileBuff(infile)) == 0) {
						fprintf(stderr, "Unexpected end of file\n");
						exit(errno | 1);
					}
					buff = infile->buffer;
				}
			}
			*num = 0;
			
			/* convert to number */
			if(mat) {
				*mat++ = strtod(strbuff, &msg);
			} else if(fmat) {
				*fmat++ = strtod(strbuff, &msg);
			} else if(smat) {
				*smat++ = dtouc(strtod(strbuff, &msg), 0.5);
			} else {
				*bmat++ = dtouc(strtod(strbuff, &msg), 0.5);
			}
			if(*msg != 0) {
				fprintf(stderr, "Malformatted entry at pos:\t(%d,%d) %s\n", dest->m, N - n, strbuff);
				exit(errno | 1);
			} else if(--avail == 0) {
				avail = buffFileBuff(infile);
				buff = infile->buffer;
			}
		}
		
		/* realloc matrix */
		if(++dest->m == dest->M) {
			Dat_realloc(dest, dest->M << 1);
			if(mat) {
				mat = dest->mat[dest->m];
			} else if(fmat) {
				fmat = dest->fmat[dest->m];
			} else if(smat) {
				smat = dest->smat[dest->m];
			} else {
				bmat = dest->bmat[dest->m];
			}
		}
	}
	
	infile->bytes = avail;
	infile->next = buff;
	
	return dest;
}

Job * loadJobs(FileBuff *infile, unsigned char sep, int col, int *N) {
	
	int avail, size, sizeJ, dim, entry, i, n, (*buffFileBuff)(FileBuff *);
	char *msg, strbuff[32];
	unsigned char c, stop, *buff, *num;
	Job *dest;
	
	/* init */
	avail = infile->bytes;
	buff = infile->next;
	buffFileBuff = infile->buffFileBuff;
	if(avail == 0) {
		if((avail = buffFileBuff(infile)) == 0) {
			return 0;
		}
		buff = infile->buffer;
	}
	
	/* skip header and get dimensions */
	do {
		dim = 1;
		while(*buff++ != '\n') {
			if(--avail == 0) {
				if((avail = buffFileBuff(infile)) == 0) {
					return 0;
				}
				buff = infile->buffer;
			}
			if(*buff == sep) {
				++dim;
			}
		}
		if(--avail == 0) {
			if((avail = buffFileBuff(infile)) == 0) {
				return 0;
			}
			buff = infile->buffer;
		}
	} while(*buff == '#');
	
	/* check dimensions */
	if(dim < col) {
		fprintf(stderr, "Invalid target column\n");
		exit(1);
	}
	
	/* allocate jobs */
	sizeJ = 1024;
	dest = job_realloc(0, 0, 0, 1024);
	n = 0; /* number of jobs */
	
	/* load rows */
	entry = 1;
	while(avail) {
		/* skip rows up to col */
		i = col;
		while(--i) {
			while((c = *buff++) != sep) {
				if(--avail == 0) {
					if((avail = buffFileBuff(infile)) == 0) {
						fprintf(stderr, "Unexpected end of file\n");
						exit(errno | 1);
					}
					buff = infile->buffer;
				}
			}
			if(--avail == 0) {
				avail = buffFileBuff(infile);
				buff = infile->buffer;
			}
		}
		
		/* read col */
		stop = dim == col ? '\n' : sep;
		num = (unsigned char *) strbuff;
		size = 32;
		while((c = *buff++) != stop) {
			*num++ = c;
			if(--size == 0) {
				fprintf(stderr, "Malformatted entry at:\t%d\n", entry);
				exit(1);
			}
			if(--avail == 0) {
				if((avail = buffFileBuff(infile)) == 0) {
					fprintf(stderr, "Unexpected end of file\n");
					exit(errno | 1);
				}
				buff = infile->buffer;
			}
		}
		*num = 0;
		if(--avail == 0) {
			avail = buffFileBuff(infile);
			buff = infile->buffer;
		}
		
		/* convert col to cluster number */
		i = strtol(strbuff, &msg, 10);
		if(*msg != 0) {
			fprintf(stderr, "Malformatted cluster at:\t%d\n", entry);
			exit(errno | 1);
		}
		
		/* check allocated jobs */
		if(sizeJ <= i) {
			size = nearestPower(i+1);
			dest = job_realloc(dest, 0, sizeJ, size);
			sizeJ = size;
		}
		
		/* update cluster entry */
		if(n <= i) {
			n = i + 1;
		}
		dest[i].size++;
		
		/* skip rest */
		if(stop != '\n') {
			while((c = *buff++) != '\n') {
				if(--avail == 0) {
					if((avail = buffFileBuff(infile)) == 0) {
						fprintf(stderr, "Unexpected end of file\n");
						exit(errno | 1);
					}
					buff = infile->buffer;
				}
			}
			if(--avail == 0) {
				avail = buffFileBuff(infile);
				buff = infile->buffer;
			}
		}
		++entry;
	}
	infile->bytes = avail;
	infile->next = buff;
	
	/* clean up jobs */
	n = cleanJobs(dest, n); /* so far n has been zero-indexed */
	if(n < sizeJ) {
		dest = job_realloc(dest, 0, sizeJ, n);
	}
	*N = n;
	
	return dest;
}

static int intCmpAscend(const void * a, const void * b) {
        return *((int*) a) - *((int*) b);
}

Job * loadMVJobs(FileBuff *infile, unsigned char sep, int col, int mv, int *MV, int *N) {
	
	int avail, size, sizeJ, dim, entry, i, j, c_i, J_i, n, nCol;
	int *MVptr, (*buffFileBuff)(FileBuff *);
	char *msg, strbuff[32];
	unsigned char c, stop, *buff, *num;
	double weight, *classes;
	Job *dest;
	
	/* init */
	avail = infile->bytes;
	buff = infile->next;
	buffFileBuff = infile->buffFileBuff;
	if(avail == 0) {
		if((avail = buffFileBuff(infile)) == 0) {
			return 0;
		}
		buff = infile->buffer;
	}
	
	/* skip header and get dimensions */
	do {
		dim = 1;
		while(*buff++ != '\n') {
			if(--avail == 0) {
				if((avail = buffFileBuff(infile)) == 0) {
					return 0;
				}
				buff = infile->buffer;
			}
			if(*buff == sep) {
				++dim;
			}
		}
		if(--avail == 0) {
			if((avail = buffFileBuff(infile)) == 0) {
				return 0;
			}
			buff = infile->buffer;
		}
	} while(*buff == '#');
	
	/* sort target columns in ascending order */
	MVptr = smalloc((mv + 1) * sizeof(int));
	*MVptr = col;
	memcpy(MVptr + 1, MV, mv * sizeof(int));
	qsort(MVptr, mv + 1, sizeof(int), intCmpAscend);
	MV = MVptr;
	
	/* check dimensions and target columns */
	if(!*MV || dim < MV[mv]) {
		fprintf(stderr, "Invalid target column\n");
		exit(1);
	} else {
		for(i = 1; i <= mv; ++i) {
			if(MV[i - 1] == MV[i]) {
				fprintf(stderr, "Invalid target columns\n");
				exit(1);
			}
		}
	}
	
	/* allocate jobs */
	sizeJ = 1024;
	dest = job_realloc(0, mv, 0, 1024);
	classes = smalloc(mv * sizeof(double));
	n = 0; /* number of jobs */
	J_i = 0;
	
	/* load rows */
	entry = 1;
	while(avail) {
		/* reset */
		c_i = 0;
		memset(classes, 0, mv * sizeof(double));
		MVptr = MV - 1;
		
		/* get columns of interest */
		stop = sep;
		nCol = mv + 2;
		j = 1;
		while(--nCol) {
			/* skip rows up to next target column */
			i = *++MVptr;
			while(i != j) {
				while((c = *buff++) != sep) {
					if(--avail == 0) {
						if((avail = buffFileBuff(infile)) == 0) {
							fprintf(stderr, "Unexpected end of file\n");
							exit(errno | 1);
						}
						buff = infile->buffer;
					}
				}
				if(--avail == 0) {
					avail = buffFileBuff(infile);
					buff = infile->buffer;
				}
				++j;
			}
			
			/* read col */
			stop = dim == i ? '\n' : sep;
			num = (unsigned char *) strbuff;
			size = 32;
			while((c = *buff++) != stop) {
				*num++ = c;
				if(--size == 0) {
					fprintf(stderr, "Malformatted entry at:\t%d\n", entry);
					exit(1);
				}
				if(--avail == 0) {
					if((avail = buffFileBuff(infile)) == 0) {
						fprintf(stderr, "Unexpected end of file\n");
						exit(errno | 1);
					}
					buff = infile->buffer;
				}
			}
			*num = 0;
			if(--avail == 0) {
				avail = buffFileBuff(infile);
				buff = infile->buffer;
			}
			++j;
			
			/* convert column to cluster or class number */
			if(i == col) {
				J_i = strtol(strbuff, &msg, 10);
				
				/* check allocated jobs */
				if(sizeJ <= J_i) {
					size = nearestPower(J_i + 1);
					dest = job_realloc(dest, mv, sizeJ, size);
					sizeJ = size;
				}
				/* update cluster entry */
				if(n <= J_i) {
					n = J_i + 1;
				}
				dest[J_i].size++;
			} else {
				weight = strtod(strbuff, &msg);
				/* update class entry */
				classes[c_i++] += weight;
			}
			if(*msg != 0) {
				fprintf(stderr, "Malformatted cluster at:\t%d\n", entry);
				exit(errno | 1);
			}
		}
		
		/* update class entries */
		for(i = 0; i < mv; ++i) {
			dest[J_i].Weights[i] += classes[i];
		}
		
		/* skip rest */
		if(stop != '\n') {
			while((c = *buff++) != '\n') {
				if(--avail == 0) {
					if((avail = buffFileBuff(infile)) == 0) {
						fprintf(stderr, "Unexpected end of file\n");
						exit(errno | 1);
					}
					buff = infile->buffer;
				}
			}
			if(--avail == 0) {
				avail = buffFileBuff(infile);
				buff = infile->buffer;
			}
		}
		++entry;
	}
	infile->bytes = avail;
	infile->next = buff;
	
	/* clean up jobs */
	n = cleanJobs(dest, n); /* so far n has been zero-indexed */
	if(n < sizeJ) {
		dest = job_realloc(dest, mv, sizeJ, n);
	}
	*N = n;
	free(MV);
	free(classes);
	
	return dest;
}

Job * loadMVEJobs(FileBuff *infile, unsigned char sep, int col, int *mv, int *N) {
	
	int avail, size, sizeJ, sizeC, dim, entry, i, j, c_i, J_i, n, nCol;
	int MV[2], *MVptr, (*buffFileBuff)(FileBuff *);
	char *msg, strbuff[32];
	unsigned char c, stop, *buff, *num;
	Job *dest;
	
	/* init */
	avail = infile->bytes;
	buff = infile->next;
	buffFileBuff = infile->buffFileBuff;
	if(avail == 0) {
		if((avail = buffFileBuff(infile)) == 0) {
			return 0;
		}
		buff = infile->buffer;
	}
	
	/* skip header and get dimensions */
	do {
		dim = 1;
		while(*buff++ != '\n') {
			if(--avail == 0) {
				if((avail = buffFileBuff(infile)) == 0) {
					return 0;
				}
				buff = infile->buffer;
			}
			if(*buff == sep) {
				++dim;
			}
		}
		if(--avail == 0) {
			if((avail = buffFileBuff(infile)) == 0) {
				return 0;
			}
			buff = infile->buffer;
		}
	} while(*buff == '#');
	
	/* sort target columns in ascending order */
	if(*mv < col) {
		MV[0] = *mv;
		MV[1] = col;
	} else {
		MV[0] = col;
		MV[1] = *mv;
	}
	
	/* check dimensions and target columns */
	if(!MV[0] || MV[0] == MV[1] || dim < MV[1]) {
		fprintf(stderr, "Invalid target column\n");
		exit(1);
	}
	
	/* allocate jobs */
	sizeJ = 1024;
	sizeC = 16;
	*mv = 0;
	dest = job_realloc(0, sizeC, 0, 1024);
	n = 0; /* number of jobs */
	J_i = 0;
	
	/* load rows */
	entry = 1;
	while(avail) {
		/* reset */
		c_i = 0;
		MVptr = MV - 1;
		
		/* get columns of interest */
		stop = sep;
		nCol = 3;
		j = 1;
		while(--nCol) {
			/* skip rows up to next target column */
			i = *++MVptr;
			while(i != j) {
				while((c = *buff++) != sep) {
					if(--avail == 0) {
						if((avail = buffFileBuff(infile)) == 0) {
							fprintf(stderr, "Unexpected end of file\n");
							exit(errno | 1);
						}
						buff = infile->buffer;
					}
				}
				if(--avail == 0) {
					avail = buffFileBuff(infile);
					buff = infile->buffer;
				}
				++j;
			}
			
			/* read col */
			stop = dim == i ? '\n' : sep;
			num = (unsigned char *) strbuff;
			size = 32;
			while((c = *buff++) != stop) {
				*num++ = c;
				if(--size == 0) {
					fprintf(stderr, "Malformatted entry at:\t%d\n", entry);
					exit(1);
				}
				if(--avail == 0) {
					if((avail = buffFileBuff(infile)) == 0) {
						fprintf(stderr, "Unexpected end of file\n");
						exit(errno | 1);
					}
					buff = infile->buffer;
				}
			}
			*num = 0;
			if(--avail == 0) {
				avail = buffFileBuff(infile);
				buff = infile->buffer;
			}
			++j;
			
			/* convert column to cluster or class number */
			if(i == col) {
				J_i = strtol(strbuff, &msg, 10);
				
				/* check allocated jobs */
				if(sizeJ <= J_i) {
					size = nearestPower(J_i + 1);
					dest = job_realloc(dest, sizeC, sizeJ, size);
					sizeJ = size;
				}
				/* update cluster entry */
				if(n <= J_i) {
					n = J_i + 1;
				}
				dest[J_i].size++;
			} else {
				c_i = strtol(strbuff, &msg, 10);
			}
			if(*msg != 0) {
				fprintf(stderr, "Malformatted cluster at:\t%d\n", entry);
				exit(errno | 1);
			}
		}
		
		/* check number of classes */
		if(sizeC <= c_i) {
			fprintf(stderr, "NONO\n");
			dest = jobWeights_realloc(dest, sizeC, sizeC << 1, sizeJ);
			sizeC <<= 1;
		}
		if(*mv <= c_i) {
			*mv = c_i + 1;
		}
		/* update class entry */
		dest[J_i].Weights[c_i] += 1;
		
		
		/* skip rest */
		if(stop != '\n') {
			while((c = *buff++) != '\n') {
				if(--avail == 0) {
					if((avail = buffFileBuff(infile)) == 0) {
						fprintf(stderr, "Unexpected end of file\n");
						exit(errno | 1);
					}
					buff = infile->buffer;
				}
			}
			if(--avail == 0) {
				avail = buffFileBuff(infile);
				buff = infile->buffer;
			}
		}
		++entry;
	}
	infile->bytes = avail;
	infile->next = buff;
	
	/* clean up jobs */
	n = cleanJobs(dest, n); /* so far n has been zero-indexed */
	if(n < sizeJ) {
		dest = job_realloc(dest, sizeC, sizeJ, n);
	}
	if(*mv < sizeC) {
		dest = jobWeights_realloc(dest, sizeC, *mv, n);
	}
	*N = n;
	
	return dest;
}
