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

Matrix * ltdMatrix_get(char *targetTemplate, char **filenames, int numFile, unsigned norm, unsigned minDepth, double (*veccmp)(short unsigned*, short unsigned*, int, int)) {
	
	/* here */
	/* get matrix of used nucs per distance */
	int i, j;
	char **filename;
	double dist, *mat;
	FileBuff *infile;
	Matrix *dest;
	MatrixCounts *mat1;
	NucCount *mat2;
	TimeStamp **targetStamp, **targetStamps;
	
	/* init */
	dest = ltdMatrix_init(numFile);
	mat1 = initMat(1048576, 128);
	mat2 = initNucCount(128);
	infile = setFileBuff(1048576);
	if(!(targetStamps = calloc(numFile, sizeof(TimeStamp *)))) {
		ERROR();
	}
	
	/* get distances */
	dest->n = numFile;
	mat = *(dest->mat);
	for(i = 0; i < numFile; ++i) {
		if((j = i) != 0) {
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
			if(FileBuffLoadMat(mat1, infile) == 0) {
				fprintf(stderr, "Malformed matrix in:\t%s\n", filenames[i]);
				exit(1);
			}
			closeFileBuff(infile);
			
			/* strip matrix for insersions */
			stripMat(mat1);
			
			/* calculate distances */
			filename = filenames;
			while(j--) {
				/* open matrix file, and find target */
				openAndDetermine(infile, *filename++);
				if(*targetStamp) {
					seekFileBiff(infile, *targetStamp);
				} else {
					while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
				}
				/* make timestamp */
				*targetStamp = timeStampFileBuff(infile, *targetStamp);
				
				/* get distance between the matrices */
				dist = cmpMats(mat1, mat2, infile, norm, minDepth, veccmp);
				if(dist < 0) {
					if(dist == -1.0) {
						fprintf(stderr, "No sufficient overlap between samples:\t%s, %s\n", filenames[i], *(filename - 1));
					} else if(dist == -2.0) {
						fprintf(stderr, "Template does not exist in file:\t%s\n", *(filename - 1));
						exit(1);
					} else {
						fprintf(stderr, "Failed to produce a distance metric between samples:\t%s, %s\n", filenames[i], *(filename - 1));
					}
				}
				*mat++ = dist;
				
				/* close mtrix file */
				closeFileBuff(infile);
				++targetStamp;
			}
		}
	}
	
	return dest;
}
