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
#include <string.h>
#include "dist.h"
#include "merge.h"
#include "tree.h"
#include "union.h"
#include "version.h"

static int helpMessage(FILE *out) {
	
	fprintf(out, "#CCPhylo enables phylogenetic analysis of samples based on overlaps between nucleotide created by e.g. KMA.\n");
	fprintf(out, "# %16s\t%-32s\n", "Options are:", "Desc:");
	fprintf(out, "# %16s\t%-32s\n", "dist", "make distance matrices");
	fprintf(out, "# %16s\t%-32s\n", "merge", "merge distance matrices");
	fprintf(out, "# %16s\t%-32s\n", "tree", "make tree(s)");
	fprintf(out, "# %16s\t%-32s\n", "union", "Find union of templates between samples");
	fprintf(out, "# %16s\t%-32s\n", "-v", "Version");
	fprintf(out, "# %16s\t%-32s\n", "-h", "Shows this helpmessage");
	return (out == stderr);
}

int main(int argc, char *argv[]) {
	
	char *arg;
	
	arg = *++argv;
	if(--argc < 1) {
		fprintf(stderr, "Too few arguments handed.\n");
		return helpMessage(stderr);
	} else if(strcmp(arg, "dist") == 0) {
		return main_dist(argc, argv);
	} else if(strcmp(arg, "tree") == 0) {
		return main_tree(argc, argv);
	} else if(strcmp(arg, "merge") == 0) {
		return main_merge(argc, argv);
	} else if(strcmp(arg, "union") == 0) {
		return main_union(argc, argv);
	} else if(strcmp(arg, "-v") == 0) {
		fprintf(stdout, "CCPhylo-%s\n", CCPHYLO_VERSION);
	} else if(strcmp(arg, "-h") == 0) {
		return helpMessage(stdout);
	} else {
		fprintf(stderr, "Unknown argument:%s\n", arg);
		return helpMessage(stderr);
	}
	
	return 0;
}
