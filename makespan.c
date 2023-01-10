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
#include "mvmakespan.h"
#include "pherror.h"
#include "tabusearch.h"
#include "mvjobs.h"
#include "mvtabusearch.h"
#include "tsv.h"

Machine * (*makespan_method)(Machine *, Job *, int, int) = &DBF;
void (*addDBEptr)(Machine **, Machine **, Job *, int, int) = &addDBE;
Machine * (*addDBFptr)(Machine *, Job *) = &addDBF;
Machine *(*FirstFitptr)(Machine *, Job *, int) = &FirstFit;
Machine *(*FirstFetptr)(Machine *, Job *) = &FirstFet;

void addDBE(Machine **Mdest, Machine **Edest, Job *J, int m, int n) {
	
	Machine *M, *E, *nextM;
	
	/* init */
	M = *Mdest;
	E = *Edest;
	
	/* put job on least loaded machine */
	M->n++;
	J->next = M->jobs;
	M->jobs = J;
	M->avail -= J->weight;
	
	/* move machine down the qeue */
	nextM = M->next;
	M->next = 0;
	if(M->n < n / m) {
		M = machinemerge(nextM, M);
	} else {
		/* remove machine from available list */
		E = machinemerge(E, M);
		M = nextM;
	}
	
	/* set pointers */
	*Mdest = M;
	*Edest = E;
}

Machine * DBE(Machine *M, Job *J, int m, int n) {
	
	Job *nextJ;
	Machine *E;
	
	/* sort both machines and jobs in descending order */
	M = machinesort(M, m);
	J = jobsort(J, n);
	
	/* init full machines */
	E = 0;
	while(J) {
		nextJ = J->next;
		
		/* put job on least loaded machine */
		if(!M) {
			/* n / m is not a Natural number */
			M = E;
			E = 0;
		}
		addDBEptr(&M, &E, J, m, n);
		
		/* move on to next job */
		J = nextJ;
	}
	
	/* merge E and M */
	M = machinemerge(M, E);
	
	return M;
}

Machine * addDBF(Machine *M, Job *J) {
	
	Machine *nextM;
	
	/* put job on least loaded machine */
	M->n++;
	J->next = M->jobs;
	M->jobs = J;
	M->avail -= J->weight;
	
	/* move machine down the qeue */
	nextM = M->next;
	M->next = 0;
	return machinemerge(nextM, M);
}

Machine * DBF(Machine *M, Job *J, int m, int n) {
	
	Job *nextJ;
	
	/* sort both machines and jobs in descending order */
	M = machinesort(M, m);
	J = jobsort(J, n);
	
	while(J) {
		nextJ = J->next;
		
		/* put job on least loaded machine */
		M = addDBFptr(M, J);
		
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
		M = FirstFitptr(M, J, m);
		
		/* move on to next job */
		J = nextJ;
	}
	
	/* un-circularize machines and return */
	nextM = M->next;
	M->next = 0;
	
	return nextM;
}

Machine * FirstFet(Machine *M, Job *J) {
	
	double weight, best;
	Machine *F, *prev, *prevF;
	
	weight = J->weight;
	best = M->avail;
	F = M;
	prev = 0;
	prevF = 0;
	while(M) {
		if(weight <= M->avail) {
			/* job fits */
			M->n++;
			J->next = M->jobs;
			M->jobs = J;
			M->avail -= weight;
			return prev;
		} else if(best < M->avail) {
			/* machine is less filled than the previous */
			best = M->avail;
			prevF = prev;
			F = M;
		}
		
		/* go to next machine */
		prev = M;
		M = M->next;
	}
	
	/* job does not fit anywhere */
	F->n++;
	J->next = F->jobs;
	F->jobs = J;
	F->avail -= weight;
	
	return prevF;
}

Machine * DFE(Machine *M, Job *J, int m, int n) {
	
	Job *nextJ;
	Machine *nextM, *E, *F;
	
	/* sort jobs in descending order */
	J = jobsort(J, n);
	
	E = 0;
	while(J) {
		nextJ = J->next;
		if(!M) {
			/* n / m is not a Natural number */
			M = E;
			E = 0;
		}
		
		/* put job on first fit (or least loaded) machine */
		F = FirstFetptr(M, J);
		
		/* remove machine from available list */
		if(F) {
			if(n / m <= F->next->n) {
				nextM = F->next;
				F->next = F->next->next;
				nextM->next = 0;
				E = machinemerge(E, nextM);
			}
		} else {
			if(n / m <= M->n) {
				nextM = M;
				M = M->next;
				nextM->next = 0;
				E = machinemerge(E, nextM);
			}
		}
		
		/* move on to next job */
		J = nextJ;
	}
	
	/* merge E and M */
	M = machinemerge(M, E);
	
	return M;
}

void print_makespan(Machine *M, FILE *out, FILE *mout) {
	
	int num, size;
	double weight;
	Job *J;
	Machine *Mptr;
	
	if(out != mout) {
		fprintf(out, "#%s\t%s\t%s\t%s\n", "Cluster", "Cluster_size", "Cluster_weight", "Partition");
		fprintf(mout, "#%s\t%s\t%s\t%s\t%s\n", "Partition", "Cluster_quantity", "Partition_size", "Partition_weight", "Partition_error");
		while(M) {
			num = M->num;
			size = 0;
			weight = 0;
			J = M->jobs;
			while(J) {
				fprintf(out, "%d\t%d\t%f\t%d\n", J->num, J->size, J->weight, num);
				size += J->size;
				weight += J->weight;
				J = J->next;
			}
			fprintf(mout, "%d\t%d\t%d\t%f\t%f\n", num, M->n, size, weight, M->avail);
			M = M->next;
		}
	} else {
		fprintf(mout, "#%s\t%s\t%s\t%s\t%s\n", "Partition", "Cluster_quantity", "Partition_size", "Partition_weight", "Partition_error");
		Mptr = M;
		while(Mptr) {
			num = Mptr->num;
			size = 0;
			weight = 0;
			J = Mptr->jobs;
			while(J) {
				size += J->size;
				weight += J->weight;
				J = J->next;
			}
			fprintf(mout, "%d\t%d\t%d\t%f\t%f\n", num, Mptr->n, size, weight, Mptr->avail);
			Mptr = Mptr->next;
		}
		
		fprintf(out, "#%s\t%s\t%s\t%s\n", "Cluster", "Cluster_size", "Cluster_weight", "Partition");
		while(M) {
			num = M->num;
			J = M->jobs;
			while(J) {
				fprintf(out, "%d\t%d\t%f\t%d\n", J->num, J->size, J->weight, num);
				J = J->next;
			}
			M = M->next;
		}
	}
}

void makespan(char *inputfilename, char *outputfilename, char *moutputfilename, int m, double *loads, int mv, int *MV, double base, unsigned char sep, int col) {
	
	int n;
	FILE *outfile, *moutfile;
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
	if(*moutputfilename == '-' && moutputfilename[1] == 0) {
		moutfile = stdout;
	} else if(strcmp(outputfilename, moutputfilename) == 0) {
		moutfile = outfile;
	} else {
		moutfile = sfopen(moutputfilename, "wb");
	}
	
	/* load jobs */
	if(!mv) {
		J = loadJobs(infile, sep, col, &n);
	} else if(MV) {
		J = loadMVJobs(infile, sep, col, mv, MV, &n);
	} else {
		J = loadMVEJobs(infile, sep, col, &mv, &n);
	}
	closeFileBuff(infile);
	destroyFileBuff(infile);
	
	/* get job weights */
	if(mv) {
		jobMVWeight(J, mv, n, base);
	} else {
		jobWeight(J, n, base);
	}
	
	/* get machines */
	if(loads) {
		M = initSkewM(m, n, mv, J, loads);
	} else {
		M = initM(m, n, mv, J);
	}
	
	/* run makespan */
	M = makespan_method(M, J, m, n);
	
	/* trade jobs */
	if(tradeM) {
		fprintf(stderr, "## Trades:\t%d\n", tradeM(M));
	}
	
	/* print results */
	print_stats(M);
	print_makespan(M, outfile, moutfile);
	if(outfile != moutfile) {
		fclose(moutfile);
	}
	fclose(outfile);
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

static int * getMVs(char *src) {
	
	int size, m, num, *dest, *dptr;
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
	dest = smalloc(++size * sizeof(int));
	dptr = dest++;
	*dptr = m;
	
	while(--size) {
		num = strtol(src, &next, 10);
		
		if(num <= 0 || (*next && *next != ',')) {
			fprintf(stderr, "Invalid multivariate cluster string:\t%s\n", next);
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
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'O', "machine_output", "Machine output file", "stdout");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'S', "separator", "Separator", "\\t");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'k', "key", "Field containing cluster number", "3");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'c', "classes", "Field(s) containing class weights", "False");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'm', "method", "Makespan initial method", "DBF");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'M', "method_help", "Help on option \"-m\"", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 't', "tabu", "Makespan tabu search method", "BB");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'T', "tabu_help", "Help on option \"-t\"", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'w', "weight", "Weighing method", "none");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'W', "weight_help", "Help on option \"-w\"", "");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'l', "loads", "Load on machines double[,double...]", "5");
	fprintf(out, "#    -%c, --%-16s\t%-32s\t%s\n", 'h', "help", "Shows this helpmessage", "");
	return (out == stderr);
}

int main_makespan(int argc, char **argv) {
	
	const char *stdstream = "-";
	int args, m, len, col, offset, mv, *MV;
	double *loads, base;
	char **Arg, *arg, *inputfilename, *outputfilename, *moutputfilename;
	char *method, *strLoads, *strMV, *trade, *weight, *errMsg, sep, opt;
	
	/* set defaults */
	m = 5;
	col = 3;
	mv = 0;
	MV = 0;
	loads = 0;
	base = 1;
	inputfilename = (char *)(stdstream);
	outputfilename = (char *)(stdstream);
	moutputfilename = (char *)(stdstream);
	method = "DBF";
	trade = "BB";
	strLoads = 0;
	strMV = 0;
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
				} else if(cmdcmp(arg, "input") == 0) {
					inputfilename = getArgDie(&Arg, &args, len + offset, "input");
				} else if(cmdcmp(arg, "output") == 0) {
					outputfilename = getArgDie(&Arg, &args, len + offset, "output");
				} else if(cmdcmp(arg, "machine_output") == 0) {
					moutputfilename = getArgDie(&Arg, &args, len + offset, "output");
				} else if(cmdcmp(arg, "separator") == 0) {
					sep = getcArgDie(&Arg, &args, len + offset, "separator");
				} else if(cmdcmp(arg, "key") == 0) {
					col = getNumArg(&Arg, &args, len + offset, "key");
				} else if(cmdcmp(arg, "classes") == 0) {
					strMV = getArgDie(&Arg, &args, len + offset, "classes");
				} else if(cmdcmp(arg, "method") == 0) {
					method = getArgDie(&Arg, &args, len + offset, "method");
				} else if(cmdcmp(arg, "method_help") == 0) {
					method = 0;
				} else if(cmdcmp(arg, "tabu") == 0) {
					trade = getArgDie(&Arg, &args, len + offset, "trade");
				} else if(cmdcmp(arg, "tabu_help") == 0) {
					trade = 0;
				} else if(cmdcmp(arg, "weight") == 0) {
					weight = getArgDie(&Arg, &args, len + offset, "weight");
				} else if(cmdcmp(arg, "weight_help") == 0) {
					weight = 0;
				} else if(cmdcmp(arg, "loads") == 0) {
					strLoads = getArgDie(&Arg, &args, len + offset, "loads");
				} else if(cmdcmp(arg, "help") == 0) {
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
					} else if(opt == 'O') {
						moutputfilename = getArgDie(&Arg, &args, len, "o");
						opt = 0;
					} else if(opt == 'S') {
						sep = getcArgDie(&Arg, &args, len, "S");
						opt = 0;
					} else if(opt == 'k') {
						col = getNumArg(&Arg, &args, len, "k");
						opt = 0;
					} else if(opt == 'c') {
						strMV = getArgDie(&Arg, &args, len, "c");
						opt = 0;
					} else if(opt == 'm') {
						method = getArgDie(&Arg, &args, len, "m");
						opt = 0;
					} else if(opt == 'M') {
						method = 0;
					} else if(opt == 't') {
						trade = getArgDie(&Arg, &args, len, "t");
						opt = 0;
					} else if(opt == 'T') {
						trade = 0;
					} else if(opt == 'w') {
						weight = getArgDie(&Arg, &args, len, "T");
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
		fprintf(stderr, "Makespan initial methods:\n");
		fprintf(stderr, "DBF:\tDecreasing Best First / Longest Processing Time (LPT)\n");
		fprintf(stderr, "DFF:\tDecreasing First Fit\n");
		fprintf(stderr, "DBE:\tDecreasing Best First with equal number of jobs\n");
		fprintf(stderr, "DFE:\tDecreasing First First with equal number of jobs\n");
		return 0;
	} else if(strcmp(method, "DBF") == 0) {
		makespan_method = &DBF;
	} else if(strcmp(method, "DFF") == 0) {
		makespan_method = &DFF;
	} else if(strcmp(method, "DBE") == 0) {
		makespan_method = &DBE;
	} else if(strcmp(method, "DFE") == 0) {
		makespan_method = &DFE;
	} else {
		invaArg("method");
	}
	
	/* get trading method */
	if(!trade) {
		fprintf(stderr, "Tabu search methods:\n");
		fprintf(stderr, "BB:\tBabettes buckets, local search + job trade\n");
		fprintf(stderr, "DBEB:\tTrades has to be with two jobs\n");
		fprintf(stderr, "None:\tNo trading\n");
		return 0;
	} else if(strcmp(trade, "BB") == 0) {
		tradeM = &tradeBB;
	} else if(strcmp(trade, "DBEB") == 0) {
		tradeM = &tradeDBEB;
	} else if(strcmp(trade, "None") == 0) {
		tradeM = 0;
	} else {
		invaArg("tabu");
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
	
	/* get multivariate classes (MV) */
	if(strMV) {
		MV = getMVs(strMV);
		mv = MV[-1];
		
		/* reset m, if balanced clusters were parsed */
		if(mv == 1) {
			mv = *MV;
			free(--MV);
			MV = 0;
		}
		if(mv < 0) {
			invaArg("classes");
		} else if(1 < mv) {
			/* set MV makespan pointers */
			addDBEptr = &addMVDBE;
			addDBFptr = &addMVDBF;
			FirstFitptr = &MVFirstFit;
			FirstFetptr = &MVFirstFet;
			negotiatePtr = &negotiateMVM;
			handoverPtr = &mvhandover;
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
		jobMVWeight = &nullMVWeight;
		base = 1;
	} else if(strncmp(weight, "log", 3) == 0) {
		jobWeight = &logWeight;
		jobMVWeight = &logMVWeight;
		base = strtod(weight + 3, &errMsg);
	} else if(strncmp(weight, "pow", 3) == 0) {
		jobWeight = &polWeight;
		jobMVWeight = &polMVWeight;
		base = strtod(weight + 3, &errMsg);
	} else if(strncmp(weight, "exp", 3) == 0) {
		jobWeight = &expWeight;
		jobMVWeight = &expMVWeight;
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
	
	/* make makespan */
	makespan(inputfilename, outputfilename, moutputfilename, m, loads, mv, MV, base, sep, col);
	
	return 0;
}
