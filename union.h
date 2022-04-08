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

#include "hashmapstr.h"

HashMapStr * unionRes(char **filenames, int numFile, char *outputfilename, double minCov, double minDepth, unsigned minLength);
int unionResPrint(char **filenames, int numFile, char *outputfilename, double minCov, double minDepth, unsigned minLength);
int unionResOrderPrint(char **filenames, int numFile, char *outputfilename, char *dbfilename, char *reffilename, double minCov, double minDepth, unsigned minLength);
int main_union(int argc, char **argv);
