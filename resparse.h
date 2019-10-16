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

#include "qseqs.h"

#ifndef RESPARSER
typedef struct resEntry ResEntry;
struct resEntry {
	Qseqs *Template;
	long unsigned Score;
	unsigned Expected;
	unsigned Template_length;
	double Template_Identity;
	double Template_Coverage;
	double Query_Identity;
	double Query_Coverage;
	double Depth;
	double q_value;
	double p_value;
};
#define RESPARSER 1
#endif

ResEntry * ResEntry_init(unsigned size);
int FileBuffValidateHeader(FileBuff *src);
int FileBuffGetEntry(FileBuff *src, ResEntry *dest);
int FileBuffGetTemplate(FileBuff *src, ResEntry *dest);
