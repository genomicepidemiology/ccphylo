/* Philip T.L.C. Clausen Jul 2022 plan@dtu.dk */

/*
 * Copyright (c) 2022, Philip Clausen, Technical University of Denmark
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
#include "cmdline.h"
#include "filebuff.h"
#include "jobs.h"
#include "machines.h"
#include "makespan.h"
#include "pherror.h"
#include "tradespan.h"
#include "tsv.h"

Machine * (*makespan_method)(Machine *M, Job *J, int m, int n) = &DBF;

Machine * DBF(Machine *M, Job *J, int m, int n) {
	
	Job *nextJ;
	Machine *nextM;
	
	/* sort both machines and jobs in descending order */
	M = machinesort(M, m);
	J = jobsort(J, n);
	
	while(J) {
		/* put job on least loaded machine */
		M->n++;
		nextJ = J->next;
		J->next = M->jobs;
		M->jobs = J;
		M->avail -= J->weight;
		
		/* move machine down the qeue */
		nextM = M->next;
		M->next = 0;
		M = machinemerge(nextM, M);
		
		/* move on to next job */
		J = nextJ;
	}
	
	return M;
}

Machine * FirstFit(Machine *M, Job *J, int m) {
	
	double weight, best;
	Machine *F;
	
	weight = J->weight;
	best = M->avail;
	F = M;
	while(m) {
		if(weight <= M->avail) {
			/* job fits */
			M->n++;
			J->next = M->jobs;
			M->jobs = J;
			M->avail -= weight;
			return M;
		} else if(best < M->avail) {
			/* machine is less filled than the previous */
			best = M->avail;
			F = M;
		}
		
		/* go to next machine */
		M = M->next;
		--m;
	}
	
	/* job does not fit anywhere */
	F->n++;
	J->next = F->jobs;
	F->jobs = J;
	F->avail -= weight;
	
	return F;
}

Machine * DFF(Machine *M, Job *J, int m, int n) {
	
	Job *nextJ;
	Machine *nextM;
	
	/* circularize machines and sort jobs in descending order */
	M[m-1].next = M;
	J = jobsort(J, n);
	
	while(J) {
		nextJ = J->next;
		
		/* put job on first fit (or least loaded) machine */
		M = FirstFit(M, J, m);
		
		/* move on to next job */
		J = nextJ;
	}
	
	/* un-circularize machines and return */
	nextM = M->next;
	M->next = 0;
	return nextM;
}

void print_makespan(Machine *M, FILE *out) {
	
	int num;
	Job *J;
	
	fprintf(out, "#%s\t%s\n", "Cluster", "Partition");
	while(M) {
		num = M->num;
		J = M->jobs;
		while(J) {
			fprintf(out, "%d\t%d\n", J->num, num);
			J = J->next;
		}
		M = M->next;
	}
}

void makespan(char *inputfilename, char *outputfilename, int m, double *loads, double base, unsigned char sep, int col) {
	
	int n;
	FILE *outfile;
	FileBuff *infile;
	Job *J;
	Machine *M;
	
	/* init file streams */
	infile = setFileBuff(1048576);
	openAndDetermine(infile, inputfilename);
	if(*outputfilename == '-' && outputfilename[1] == 0) {
		outfile = stdout;
	} else {
		outfile = sfopen(outputfilename, "wb");
	}
	
	/* load jobs */
	J = loadJobs(infile, sep, col, &n);
	closeFileBuff(infile);
	destroyFileBuff(infile);
	
	/* get job weights */
	jobWeight(J, n, base);
	
	/* rm invalid jobs (weight <= 0) */
	//n = cleanJobs(J, n);
	
	/* get machines */
	if(loads) {
		M = initSkewM(m, n, J, loads);
	} else {
		M = initM(m, n, J);
	}
	
	/* run makespan */
	M = makespan_method(M, J, m, n);
	
	/* trade jobs */
	fprintf(stderr, "Number of trades:\t%d\n", tradeM(M));
	
	/* print results */
	print_makespan(M, outfile);
	fclose(outfile);
	fprintf(stderr, "MSE:\t%f\n", machineMSE(M));
}

static double * getLoads(char *src) {
	
	int size, m;
	double num, *dest, *dptr;
	char *ptr, *next, c;
	
	/* get size */
	size = 1;
	ptr = src;
	while((c = *++ptr)) {
		if(c == ',') {
			++size;
		}
	}
	
	m = size;
	dest = smalloc(++size * sizeof(double));
	dptr = dest++;
	*dptr = m;
	
	while(--size) {
		num = strtod(src, &next);
		
		if(num <= 0 || (*next && *next != ',')) {
			fprintf(stderr, "Invalid load string:\t%s\n", next);
			exit(1);
		} else {
			*++dptr = num;
		}
		src = next + 1;
	}
	
	return dest;
}

static int helpMessage(FILE *out) {
	
	fprintf(out, "#CCPhylo make a DBSCAN given a set of phylip distance matrices.\n");
	fprintf(out, "#   %-24s\t%-32s\t%s\n", "Options are:", "Desc:", "Default:");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'i', "input", "Input file", "stdin");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'o', "output", "Output file", "stdout");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'S', "separator", "Separator", "\\t");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'k', "key", "Field containing cluster number", "3");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'm', "method", "Makespan method", "DBF");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'M', "method_help", "Help on option \"-m\"", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'w', "weight", "Weighing method", "none");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'W', "weight_help", "Help on option \"-w\"", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'l', "loads", "Load on machines double[,double...]", "5");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'h', "help", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_makespan(int argc, char **argv) {
	
	const char *stdstream = "-";
	int args, m, len, col, offset;
	double *loads, base;
	char **Arg, *arg, *inputfilename, *outputfilename, *method, *strLoads;
	char *weight, *errMsg, sep, opt;
	
	/* set defaults */
	m = 5;
	col = 3;
	loads = 0;
	base = 1;
	inputfilename = (char *)(stdstream);
	outputfilename = (char *)(stdstream);
	method = "DBF";
	strLoads = 0;
	weight = "none";
	errMsg = 0;
	sep = '\t';
	
	/* parse cmd-line */
	args = argc - 1;
	Arg = argv;
	if(args && **++Arg == '-') {
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
				offset = 2 + (arg[len] ? 1 : 0);
				
				/* long option */
				if(*arg == 0) {
					/* terminate cmd-line */
					++Arg;
				} else if(strncmp(arg, "input", len) == 0) {
					inputfilename = getArgDie(&Arg, &args, len + offset, "input");
				} else if(strncmp(arg, "output", len) == 0) {
					outputfilename = getArgDie(&Arg, &args, len + offset, "output");
				} else if(strncmp(arg, "separator", len) == 0) {
					sep = getcArgDie(&Arg, &args, len + offset, "separator");
				} else if(strncmp(arg, "key", len) == 0) {
					col = getNumArg(&Arg, &args, len + offset, "key");
				} else if(strncmp(arg, "method", len) == 0) {
					method = getArgDie(&Arg, &args, len + offset, "method");
				} else if(strncmp(arg, "method_help", len) == 0) {
					method = 0;
				} else if(strncmp(arg, "weight", len) == 0) {
					weight = getArgDie(&Arg, &args, len + offset, "weight");
				} else if(strncmp(arg, "weight_help", len) == 0) {
					weight = 0;
				} else if(strncmp(arg, "loads", len) == 0) {
					strLoads = getArgDie(&Arg, &args, len + offset, "loads");
				} else if(strncmp(arg, "help", len) == 0) {
					return helpMessage(stdout);
				} else {
					unknArg(arg - 2);
				}
			} else {
				/* multiple option */
				len = 1;
				opt = *arg;
				while(opt && (opt = *arg++)) {
					++len;
					if(opt == 'i') {
						inputfilename = getArgDie(&Arg, &args, len, "i");
						opt = 0;
					} else if(opt == 'o') {
						outputfilename = getArgDie(&Arg, &args, len, "o");
						opt = 0;
					} else if(opt == 'S') {
						sep = getcArgDie(&Arg, &args, len, "S");
						opt = 0;
					} else if(opt == 'k') {
						col = getNumArg(&Arg, &args, len, "k");
						opt = 0;
					} else if(opt == 'm') {
						method = getArgDie(&Arg, &args, len, "m");
						opt = 0;
					} else if(opt == 'M') {
						method = 0;
					} else if(opt == 'w') {
						weight = getArgDie(&Arg, &args, len, "w");
						opt = 0;
					} else if(opt == 'W') {
						weight = 0;
					} else if(opt == 'l') {
						strLoads = getArgDie(&Arg, &args, len, "l");
						opt = 0;
					} else if(opt == 'h') {
						return helpMessage(stdout);
					} else {
						*arg = 0;
						unknArg(arg - 1);
					}
				}
			}
		} else {
			/* terminate cmd-line */
			--arg;
			++args;
			len = 0;
		}
		--args;
	}
	
	/* non-options */
	if(args) {
		inputfilename = *Arg;
		if(--args) {
			nonOptError();
		}
	}
	
	/* get method */
	if(!method) {
		fprintf(stderr, "Makespan methods:\n");
		fprintf(stderr, "DBF:\tDecreasing Best First\n");
		fprintf(stderr, "DFF:\tDecreasing First Fit\n");
		return 0;
	} else if(strcmp(method, "DBF") == 0) {
		makespan_method = &DBF;
	} else if(strcmp(method, "DFF") == 0) {
		makespan_method = &DFF;
	} else {
		invaArg("method");
	}
	
	/* get loads */
	if(strLoads) {
		loads = getLoads(strLoads);
		m = loads[-1];
		
		/* reset m, if balanced load were parsed */
		if(m == 1) {
			m = *loads;
			free(--loads);
			loads = 0;
		}
		if(m <= 0) {
			invaArg("loads");
		}
	}
	
	/* get weight */
	if(!weight) {
		fprintf(stderr, "Weight methods:\n");
		fprintf(stderr, "none:\tDo not weigh clusters\n");
		fprintf(stderr, "logX:\tWeigh one plus logarithmicly with base X\n");
		fprintf(stderr, "powX:\tWeigh polynomial with exponent X\n");
		fprintf(stderr, "expX:\tWeigh exponential with exponential base X\n");
		return 0;
	} else if(strcmp(weight, "none") == 0) {
		jobWeight = &nullWeight;
		base = 1;
	} else if(strncmp(weight, "log", 3) == 0) {
		jobWeight = &logWeight;
		base = strtod(weight + 3, &errMsg);
	} else if(strncmp(weight, "pow", 3) == 0) {
		jobWeight = &polWeight;
		base = strtod(weight + 3, &errMsg);
	} else if(strncmp(weight, "exp", 3) == 0) {
		jobWeight = &expWeight;
		base = strtod(weight + 3, &errMsg);
	} else {
		invaArg("weight");
	}
	if(errMsg && *errMsg) {
		if(strcmp(errMsg, "e") == 0) {
			base = 2.71828182845904523536028747135266;
		} else {
			invaArg("weight");
		}
	}
	
	/* make tree */
	makespan(inputfilename, outputfilename, m, loads, base, sep, col);
	
	return 0;
}
