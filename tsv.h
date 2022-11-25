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

#include "dat.h"
#include "filebuff.h"
#include "jobs.h"

Dat * loadTsv(FileBuff *infile, unsigned char sep);
Job * loadJobs(FileBuff *infile, unsigned char sep, int col, int *N);
Job * loadMVJobs(FileBuff *infile, unsigned char sep, int col, int mv, int *MV, int *N);
Job * loadMVEJobs(FileBuff *infile, unsigned char sep, int col, int *mv, int *N);
