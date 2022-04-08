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
#include "pherror.h"
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
