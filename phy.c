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

#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include "filebuff.h"
#include "matrix.h"
#include "pherror.h"
#include "phy.h"
#include "qseqs.h"

static char * stripDir(char *str) {
	
	char *ptr;
	
	/* strip from leading directories */
	ptr = str - 1;
	while(*++ptr) {
		if(*ptr == '/') {
			str = ptr + 1;
		}
	}
	
	return str;
}

void printphy(FILE *outfile, Matrix *src, char **names, unsigned format) {
	
	int i, j;
	double *ptr;
	
	/* here */
	Qseqs **Names = (Qseqs **) names;
	
	fprintf(outfile, "%10d\n", src->n);
	ptr = *(src->mat);
	for(i = 0; i < src->n; ++i) {
		if(format & 1) {
			fprintf(outfile, "%s", Names[i]->seq);
		} else {
			fprintf(outfile, "%10s", stripDir(names[i]));
		}
		
		if((j = i) != 0) {
			while(j--) {
				fprintf(outfile, "\t%.4f", *ptr++);
			}
		}
		fprintf(outfile, "\n");
	}
}

Qseqs ** loadPhy(Matrix *src, Qseqs **names, FileBuff *infile) {
	
	int i, j, n, avail, size, (*buffFileBuff)(FileBuff *);
	char *msg, strbuff[32];
	unsigned char c, stop, *buff, *seq;
	double *mat;
	Qseqs *name;
	
	/* init */
	avail = infile->bytes;
	buff = infile->next;
	buffFileBuff = infile->buffFileBuff;
	if(avail == 0) {
		if((avail = buffFileBuff(infile)) == 0) {
			src->n = 0;
			return names;
		}
		buff = infile->buffer;
	}
	
	/* get matrix size */
	n = 0;
	while((c = *buff++) != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(infile)) == 0) {
				src->n = 0;
				return names;
			}
			buff = infile->buffer;
		}
		if('0' <= c && c <= '9') {
			n = 10 * n + (c - '0');
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(infile)) == 0) {
			src->n = 0;
			return names;
		}
		buff = infile->buffer;
	}
	
	/* alloc */
	if(!names) {
		if(src->size < n) {
			ltdMatrix_realloc(src, n);
		}
		names = smalloc(n * sizeof(Qseqs *));
		i = n;
		while(i--) {
			names[i] = setQseqs(32);
		}
	} else if(src->size < n) {
		i = src->size - 1;
		ltdMatrix_realloc(src, n);
		if(!(names = realloc(names, n * sizeof(Qseqs *)))) {
			ERROR();
		}
		while(++i < n) {
			names[i] = setQseqs(32);
		}
	}
	
	/* validate size */
	if(!(src->n = n)) {
		return names;
	}
	
	/* load rows */
	mat = *(src->mat);
	for(i = 0; i < n; ++i) {
		/* get name */
		name = names[i];
		seq = name->seq;
		size = name->size;
		while((c = (*seq++ = *buff++)) != '\t' && c != '\n') {
			if(--avail == 0) {
				if((avail = buffFileBuff(infile)) == 0) {
					fprintf(stderr, "Malformatted phylip file, name on row: %d\n", ++i);
					errno |= 1;
					return 0;
				}
				buff = infile->buffer;
			}
			if(--size == 0) {
				size = name->size;
				name->size <<= 1;
				name->seq = realloc(name->seq, name->size);
				if(!name->seq) {
					ERROR();
				}
				seq = name->seq + size;
			}
		}
		/* chomp seq */
		*seq = 0;
		while(isspace(*--seq)) {
			*seq = 0;
			++size;
		}
		name->len = name->size - size;
		if(--avail == 0) {
			if((avail = buffFileBuff(infile)) == 0) {
				fprintf(stderr, "Malformatted phylip file, name on row: %d\n", ++i);
				errno |= 1;
				return 0;
			}
			buff = infile->buffer;
		}
		
		/* get distances */
		j = i;
		while(j--) {
			stop = j != 0 ? '\t' : '\n';
			seq = (unsigned char *) strbuff;
			*seq = 0;
			while((c = *buff++) != stop && c != '\t') {
				*seq++ = c;
				if(--avail == 0) {
					if((avail = buffFileBuff(infile)) == 0) {
						++i;
						fprintf(stderr, "Malformatted phylip file, distance pos:\t(%d,%d)\n", i, i - j);
						errno |= 1;
						return 0;
					}
					buff = infile->buffer;
				}
			}
			*seq = 0;
			
			*mat++ = strtod(strbuff, &msg);
			if(*msg != 0) {
				++i;
				fprintf(stderr, "Malformatted distance as pos:\t(%d,%d)\n", i, i - j);
				exit(errno | 1);
			} else if(--avail == 0 && (stop != '\n' || i != n - 1)) {
				if((avail = buffFileBuff(infile)) == 0) {
					return 0;
				}
				buff = infile->buffer;
			}
		}
		
		while(c != '\n') {
			c = *buff++;
			if(--avail == 0) {
				if((avail = buffFileBuff(infile)) == 0 && i != n - 1) {
					fprintf(stderr, "Malformatted phylip file, missing newline at row:\t%d\n", ++i);
					errno |= 1;
					return 0;
				}
				buff = infile->buffer;
			}
		}
	}
	
	infile->bytes = avail;
	infile->next = buff;
	
	return names;
}
