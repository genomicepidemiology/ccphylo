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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fbseek.h"
#include "filebuff.h"
#include "matcmp.h"
#include "matparse.h"
#include "matrix.h"
#include "pherror.h"
#define strcmp2(str1, str2) (*(str1) == *(str2) && strcmp(str1 + 1, str2 + 1) == 0)

void ltdMatrix_get(Matrix *dest, Matrix *nDest, MatrixCounts *mat1, NucCount *mat2, FileBuff *infile, TimeStamp **targetStamps, unsigned char *include, char *targetTemplate, char **filenames, int numFile, unsigned norm, unsigned minDepth, unsigned minLength, double minCov, double (*veccmp)(short unsigned*, short unsigned*, int, int)) {
	
	int i, j, n;
	char **filename;
	double dist, *mat, *nMat;
	TimeStamp **targetStamp;
	
	/* get distances */
	dest->n = 0;
	mat = *(dest->mat);
	nMat = *(nDest->mat);
	for(i = 0; i < numFile; ++i) {
		if((j = i) != 0) {
			if(include[i]) {
				/* open matrix file, and find target */
				openAndDetermine(infile, filenames[i]);
				if(*(targetStamp = targetStamps)) {
					seekFileBiff(infile, *targetStamp);
				} else {
					while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
				}
				/* make timestamp */
				*targetStamp = timeStampFileBuff(infile, *targetStamp);
				
				/* initialize mat1 */
				setMatName(mat1, mat2);
				if(FileBuffLoadMat(mat1, infile, minDepth) == 0) {
					fprintf(stderr, "Malformed matrix in:\t%s\n", filenames[i]);
					exit(1);
				}
				closeFileBuff(infile);
				
				/* validate matrix */
				if(mat1->nNucs < minLength || mat1->nNucs < minCov * mat1->len) {
					fprintf(stderr, "Template did exceed threshold for inclusion:\t%s\n", filenames[i]);
					include[i] = 0;
				}
				
				/* strip matrix for insersions */
				stripMat(mat1);
			}
			
			if(include[i]) {
				/* calculate distances */
				filename = filenames;
				n = 0;
				while(j--) {
					if(include[n]) {
						/* open matrix file, and find target */
						openAndDetermine(infile, *filename);
						if(*targetStamp) {
							seekFileBiff(infile, *targetStamp);
						} else {
							while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
						}
						/* make timestamp */
						*targetStamp = timeStampFileBuff(infile, *targetStamp);
						
						/* get distance between the matrices */
						dist = cmpMats(mat1, mat2, infile, norm, minDepth, minLength, minCov, veccmp);
						if(dist < 0) {
							if(dist == -1.0) {
								fprintf(stderr, "No sufficient overlap between samples:\t%s, %s\n", filenames[i], *filename);
							} else if(dist == -2.0) {
								fprintf(stderr, "Template did exceed threshold for inclusion:\t%s\n", *filename);
							} else {
								fprintf(stderr, "Failed to produce a distance metric between samples:\t%s, %s\n", filenames[i], *filename);
							}
						}
						
						if(-1.0 <= dist) {
							*mat++ = dist;
							*nMat++ = mat2->total;
						} else {
							include[n] = 0;
						}
						
						/* close mtrix file */
						closeFileBuff(infile);
						++targetStamp;
					}
					++filename;
					++n;
				}
			}
		}
	}
	
	/* get number of include templates */
	++numFile;
	--include;
	n = 0;
	while(--numFile) {
		if(*++include) {
			++n;
		}
	}
	dest->n = n;
	nDest->n = n;
	
}

int ltdRow_get(double *D, double *N, MatrixCounts *mat1, NucCount *mat2, FileBuff *infile, char *targetTemplate, char *addfilename, Qseqs **filenames, int n, unsigned norm, unsigned minDepth, unsigned minLength, double minCov, double (*veccmp)(short unsigned*, short unsigned*, int, int)) {
	
	char *filename;
	double dist, *mat, *nMat;
	
	/* load new sample matrix into memory */
	/* open matrix file, and find target */
	openAndDetermine(infile, addfilename);
	while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
	
	/* initialize mat1 */
	setMatName(mat1, mat2);
	if(FileBuffLoadMat(mat1, infile, minDepth) == 0) {
		fprintf(stderr, "Malformed matrix in:\t%s\n", addfilename);
		exit(1);
	}
	closeFileBuff(infile);
	
	/* validate matrix */
	if(mat1->nNucs < minLength || mat1->nNucs < minCov * mat1->len) {
		fprintf(stderr, "Template did exceed threshold for inclusion:\t%s\n", addfilename);
		return 1;
	}
	
	/* strip matrix for insersions */
	stripMat(mat1);
	
	/* calculate distances */
	mat = D;
	nMat = N;
	--filenames;
	++n;
	while(--n) {
		filename = (char *) (*++filenames)->seq;
		/* open matrix file, and find target */
		openAndDetermine(infile, filename);
		while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
		
		/* get distance between the matrices */
		dist = cmpMats(mat1, mat2, infile, norm, minDepth, minLength, minCov, veccmp);
		if(dist < 0) {
			if(dist == -1.0) {
				fprintf(stderr, "No sufficient overlap between samples:\t%s, %s\n", filename, addfilename);
			} else if(dist == -2.0) {
				fprintf(stderr, "Template did exceed threshold for inclusion:\t%s\n", filename);
				return 1;
			} else {
				fprintf(stderr, "Failed to produce a distance metric between samples:\t%s, %s\n", filename, addfilename);
				return 1;
			}
		}
		
		*mat++ = dist;
		*nMat++ = mat2->total;
		
		/* close mtrix file */
		closeFileBuff(infile);
	}
	
	return 0;
}
