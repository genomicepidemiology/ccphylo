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
#include <string.h>
#include "filebuff.h"
#include "matrix.h"
#include "pherror.h"
#include "phy.h"
#include "qseqs.h"

char * stripDir(char *str) {
	
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

void printphy(FILE *outfile, Matrix *src, char **names, unsigned char *include, char *comment, unsigned format) {
	
	int i, j, jStart;
	double d, *ptr;
	
	/* printf comment */
	if(format & 4) {
		fprintf(outfile, "#%s\n", comment);
	}
	fprintf(outfile, "%10d\n", src->n);
	ptr = *(src->mat);
	jStart = 0;
	for(i = 0; jStart != src->n; ++i) {
		if(include == 0 || include[i]) {
			if(format & 1) {
				fprintf(outfile, "%s", stripDir(names[i]));
			} else {
				fprintf(outfile, "%-10.10s", stripDir(names[i]));
			}
			
			j = jStart++;
			while(j--) {
				d = *ptr++;
				if(d == (int) d) {
					fprintf(outfile, "\t%d", (int) d);
				} else {
					fprintf(outfile, "\t%.4f", d);
				}
			}
			fprintf(outfile, "\n");
		}
	}
}

void printphyUpdate(FILE *outfile, int n, char *name, double *D, unsigned format) {
	
	int c;
	double d;
	
	/* skip comment */
	if((c = getc(outfile)) == '#') {
		while((c = getc(outfile)) != '\n' && c != EOF);
	} else {
		ungetc(c, outfile);
	}
	/* print new size */
	fseek(outfile, 0, SEEK_SET);
	fprintf(outfile, "%10d", n);
	fflush(outfile);
	
	/* print new row */
	fseek(outfile, 0, SEEK_END);
	if(format & 1) {
		fprintf(outfile, "%s", stripDir(name));
	} else {
		fprintf(outfile, "%-10.10s", stripDir(name));
	}
	--D;
	while(--n) {
		d = *++D;
		if(d == (int) d) {
			fprintf(outfile, "\t%d", (int) d);
		} else {
			fprintf(outfile, "\t%.4f", d);
		}
	}
	fprintf(outfile, "\n");
	fflush(outfile);
}

Qseqs ** loadPhy(Matrix *src, Qseqs **names, Qseqs *header, FileBuff *infile) {
	
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
	
	/* get comment */
	if(*buff == '#') {
		if(header) {
			size = header->size;
			seq = header->seq - 1;
			++buff;
			if(--avail == 0) {
				if((avail = buffFileBuff(infile)) == 0) {
					src->n = 0;
					return names;
				}
				buff = infile->buffer;
			}
			while((*++seq = *buff++) != '\n') {
				if(--avail == 0) {
					if((avail = buffFileBuff(infile)) == 0) {
						src->n = 0;
						return names;
					}
					buff = infile->buffer;
				}
				if(--size == 0) {
					size = header->size;
					header->size <<= 1;
					if(!(header->seq = realloc(header->seq, header->size))) {
						ERROR();
					}
					seq = header->seq + size;
				}
			}
			*seq = 0;
			header->len = header->size - size;
			
			if(--avail == 0) {
				if((avail = buffFileBuff(infile)) == 0) {
					src->n = 0;
					return names;
				}
				buff = infile->buffer;
			}
		} else {
			while(*buff++ != '\n') {
				if(--avail == 0) {
					if((avail = buffFileBuff(infile)) == 0) {
						src->n = 0;
						return names;
					}
					buff = infile->buffer;
				}
			}
			if(--avail == 0) {
				if((avail = buffFileBuff(infile)) == 0) {
					src->n = 0;
					return names;
				}
				buff = infile->buffer;
			}
		}
	} else if(header) {
		header->len = 0;
		*(header->seq) = 0;
		
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
		names = smalloc(src->size * sizeof(Qseqs *));
		i = src->size;
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
					src->n = 0;
					return names;
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
		name->len = name->size - size + 1;
		if(--avail == 0) {
			if((avail = buffFileBuff(infile)) == 0) {
				fprintf(stderr, "Malformatted phylip file, name on row: %d\n", ++i);
				errno |= 1;
				src->n = 0;
				return names;
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
						src->n = 0;
						return names;
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
					src->n = 0;
					return names;
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
					src->n = 0;
					return names;
				}
				buff = infile->buffer;
			}
		}
	}
	
	infile->bytes = avail;
	infile->next = buff;
	
	return names;
}

int getSizePhy(FileBuff *infile) {
	
	int n, avail, (*buffFileBuff)(FileBuff *);
	unsigned char c, *buff;
	
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
	
	/* skip comment */
	if(*buff == '#') {
		while(*buff++ != '\n') {
			if(--avail == 0) {
				if((avail = buffFileBuff(infile)) == 0) {
					return 0;
				}
				buff = infile->buffer;
			}
		}
		if(--avail == 0) {
			if((avail = buffFileBuff(infile)) == 0) {
				return 0;
			}
			buff = infile->buffer;
		}
	}
	
	/* get matrix size */
	n = 0;
	while((c = *buff++) != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(infile)) == 0) {
				return 0;
			}
			buff = infile->buffer;
		}
		if('0' <= c && c <= '9') {
			n = 10 * n + (c - '0');
		}
	}
	
	/* set infile */
	infile->bytes = avail - 1;
	infile->next = buff;
	
	return n;
}

Qseqs ** getFilenamesPhy(char *path, int n, FileBuff *infile) {
	
	int i, avail, size, len, (*buffFileBuff)(FileBuff *);
	unsigned char c, *buff, *seq;
	Qseqs *name, **names;
	
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
	
	/* alloc and init names */
	len = strlen(path);
	names = smalloc(n * sizeof(Qseqs *));
	i = n;
	while(i--) {
		names[i] = setQseqs(len + 32);
		names[i]->len = len;
		memcpy(names[i]->seq, path, len);
	}
	
	/* load names */
	for(i = 0; i < n; ++i) {
		/* get name */
		name = names[i];
		seq = name->seq + len;
		size = name->size - len;
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
				seq = name->seq + len + size;
			}
		}
		/* chomp seq */
		*seq = 0;
		while(isspace(*--seq)) {
			*seq = 0;
			++size;
		}
		name->len = name->size - size + 1;
		if(--avail == 0) {
			if((avail = buffFileBuff(infile)) == 0) {
				fprintf(stderr, "Malformatted phylip file, name on row: %d\n", ++i);
				errno |= 1;
				return 0;
			}
			buff = infile->buffer;
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
