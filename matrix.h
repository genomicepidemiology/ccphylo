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

#ifndef MATRIX
typedef struct matrix Matrix;
struct matrix {
	int n;
	int size;
	double **mat;
};
#define MATRIX 1
#endif

Matrix * ltdMatrix_init(unsigned size);
void ltdMatrix_realloc(Matrix *src, unsigned size);
void Matrix_destroy(Matrix *src);
void ltdMatrix_popArrange(Matrix *mat, unsigned pos);
