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
#include "cmdline.h"
#include "dbscan.h"
#include "dist.h"
#include "fullphy.h"
#include "makespan.h"
#include "merge.h"
#include "nwck2phy.h"
#include "phycmp.h"
#include "rarify.h"
#include "tree.h"
#include "trim.h"
#include "tsv2phy.h"
#include "union.h"
#include "version.h"

static int helpMessage(FILE *out) {
	
	/* All options:
	a	add
	A	fragment_amount
	b	byte_precision
	B	database
	c	citation
	C	min_cov
	d	distance
	D	distance_help
	e	max_distance
	E	min_depth
	f	flag
	F	flag_help
	g	free
	h	help
	H	mmap
	i	input
	l	significance_lvl
	L	min_len
	m	method
	M	method_help
	n	nucleotide_numbers
	N	min_neighbors
	o	output
	O	nucletides_included
	P	proximity
	p	float_precision
	r	reference
	R	rarification_factor
	s	short_precision 
	S	separator
	t	threads
	T	tmp
	v	version
	V	nucleotide_variations
	w	nucleotides_weights
	W	normalization_weight
	y	methylation_motifs
	*/
	
	fprintf(out, "# CCPhylo enables phylogenetic analysis of samples based on overlaps between nucleotide created by e.g. KMA. Input file(s) may be given as non-option arguments succeding all options\n");
	fprintf(out, "#   %-24s\t%-32s\n", "Options are:", "Desc:");
	fprintf(out, "#    %-23s\t%-32s\n", "dist", "make distance matrices");
	fprintf(out, "#    %-23s\t%-32s\n", "tree", "make tree(s)");
	fprintf(out, "#    %-23s\t%-32s\n", "dbscan", "make DBSCAN(s)");
	fprintf(out, "#    %-23s\t%-32s\n", "union", "Find union of templates between samples");
	fprintf(out, "#    %-23s\t%-32s\n", "merge", "merge distance matrices");
	fprintf(out, "#    %-23s\t%-32s\n", "nwck2phy", "Convert newick file to phylip distance file");
	fprintf(out, "#    %-23s\t%-32s\n", "tsv2phy", "Convert tsv file to phylip distance file");
	fprintf(out, "#    %-23s\t%-32s\n", "rarify", "Rarify a KMA matrix");
	fprintf(out, "#    %-23s\t%-32s\n", "trim", "Trim multiple alignments");
	fprintf(out, "#    %-23s\t%-32s\n", "phycmp", "Compare two phylip matrices");
	fprintf(out, "#    %-23s\t%-32s\n", "fullphy", "Convert ltd phy to full phy");
	fprintf(out, "#    %-23s\t%-32s\n", "makespan", "make Makespan");
	fprintf(out, "#    -%c, --%-17s\t%-32s\n", 'v', "version", "Version");
	fprintf(out, "#    -%c, --%-17s\t%-32s\n", 'c', "citation", "Citation");
	fprintf(out, "#    -%c, --%-17s\t%-32s\n", 'h', "help", "Shows this helpmessage");
	return (out == stderr);
}

int main(int argc, char **argv) {
	
	int args, len, flag;//, offset;
	char **Arg, *arg, opt;
	
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
	} else if(strcmp(arg, "nwck2phy") == 0) {
		return main_nwck2phy(argc, argv);
	} else if(strcmp(arg, "tsv2phy") == 0) {
		return main_tsv2phy(argc, argv);
	} else if(strcmp(arg, "rarify") == 0 || strcmp(arg, "rarefy") == 0) {
		return main_rarify(argc, argv);
	} else if(strcmp(arg, "trim") == 0) {
		return main_trim(argc, argv);
	} else if(strcmp(arg, "dbscan") == 0) {
		return main_dbscan(argc, argv);
	} else if(strcmp(arg, "phycmp") == 0) {
		return main_phycmp(argc, argv);
	} else if(strcmp(arg, "fullphy") == 0) {
		return main_fullphy(argc, argv);
	} else if(strcmp(arg, "makespan") == 0) {
		return main_makespan(argc, argv);
	} else if(*arg == '-') {
		/* parse options */
		flag = 0;
		args = argc;
		Arg = argv;
		if(args && **Arg == '-') {
			len = 1;
			--Arg;
		} else {
			len = 0;
		}
		while(args && len) {
			arg = *++Arg;
			if(*arg++ == '-') {
				if(*arg == '-') {
					/* check if argument is included */
					len = getOptArg(++arg);
					//offset = 2 + (arg[len] ? 1 : 0);
					
					/* long option */
					if(*arg == 0) {
						/* terminate cmd-line */
						++Arg;
					} else if(cmdcmp(arg, "version") == 0) {
						flag |= 1;
					} else if(cmdcmp(arg, "citation") == 0) {
						flag |= 2;
					} else if(cmdcmp(arg, "help") == 0) {
						flag |= 4;
					} else {
						unknArg(arg - 2);
					}
				} else {
					/* multiple option */
					len = 1;
					opt = *arg;
					while(opt && (opt = *arg++)) {
						++len;
						if(opt == 'v') {
							flag |= 1;
						} else if(opt == 'c') {
							flag |= 2;
						} else if(opt == 'h') {
							flag |= 4;
						} else {
							*arg = 0;
							unknArg(arg - 1);
						}
					}
				}
			} else {
				unknArg(arg - 1);
			}
			--args;
		}
		
		/* non-options */
		if(args) {
			nonOptError();
		}
		
		/* print specified output */
		if(flag & 1) {
			fprintf(stdout, "CCPhylo-%s\n", CCPHYLO_VERSION);
		}
		if(flag & 2) {
			fprintf(stdout, "1. Philip T.L.C. Clausen, \"Scaling neighbor joining to one million taxa with dynamic and heuristic neighbor joining\", Bioinformatics, 2023, https://doi.org/10.1093/bioinformatics/btac774.\n");
			fprintf(stdout, "2. Malte B. Hallgren, Soeren Overballe-Petersen, Ole Lund, Henrik Hasman, Philip T.L.C. Clausen, \"MINTyper: an outbreak-detection method for accurate and rapid SNP typing of clonal clusters with noisy long reads\", Biology Methods & Protocols, 2021, https://doi.org/10.1093/biomethods/bpab008.\n");
		}
		if(flag & 4) {
			return helpMessage(stdout);
		}
	} else {
		fprintf(stderr, "Unknown argument:%s\n", arg);
		return helpMessage(stderr);
	}
	
	return 0;
}
