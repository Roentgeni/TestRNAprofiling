#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hashtbl.h"
#include "Set.h"
#include "Options.h"
#include "memoryDFS.h"
#include "boltzmann_main.h"
#include <math.h>
#include <time.h>
// #include <sys/time.h>

using namespace std;
int main0(int argc, char *argv[],int flag, FILE* data,FILE* re,int writedata);
int copy(char* to, char* from, int start, int end );

int copy(char* to, char* from, int start, int end ){
  int t=0;
  if(end>=strlen(from)){
    end=strlen(from);
  }
  for(int i=start;i<=end;i++){
    to[t]=from[i];
    if(from[i]==' ')break;
    t++;
  }
  to[t]='\0';
  return 0;
}

int main(int argc, char *argv[]){
  
  int n=atoi(argv[argc-1]);
  FILE* r;//source file
  FILE* result;
  FILE* not0;
  FILE* data;
  int writedata=0;
  clock_t time;
  clock_t totaltime;
  totaltime=clock();
  r=fopen(argv[argc-2],"r");
  result=fopen("temp.txt","w+");
  not0=fopen("delete.txt","w+");
  data=fopen("data.txt","w+");
    if(r==NULL){
      printf("cannot open file r %s\n",argv[argc-2]);
      return 0;
    }
    if(result==NULL){
      printf("%s\n","cannot open file result");
      return 0;
    }
    if(not0==NULL){
      printf("%s\n","cannot open file not0(delete.txt)");
      return 0;
    }
    if(data==NULL){
      printf("%s\n","cannot open file data");
      return 0;
    }
    char name[256];
    fscanf(r, "%[^\n] ", name);
    while(1) {writedata++;
      time=clock();
        result=fopen("temp.txt","a+");
		if(name[0]=='>'){
		  //create a new fasta file for the RNA
			char filename[25];
			copy(filename,name,1,21);
			FILE* w;
			w=fopen(filename,"w");
			if(w==NULL){
			  printf("cannot open file w\n");
			  return 0;
			}
			fprintf(w,"%s\n",name);
			//fwrite(name,1,strlen(name),w);fwrite("\n",1,1,w);
			while(1) {
			  char buf[100];
			  fscanf(r, "%[^\n] ", buf);
			  if(buf[0]=='>'){
			    memset(name, '\0', sizeof(name));
			    strcpy(name,buf);
			    break;
			  }
			  // fwrite(buf,1,strlen(buf),w);fwrite("\n",1,1,w);
			  fprintf(w,"%s\n",buf);
			  if(feof(r)){
			    break;
			  }
			  
			}
			fclose(w);//finish writing the fasta file
		        
			char* argvs[2];
			argvs[1]=filename;
			int x[n];
			int i;
			float average, variance, std_deviation, sum = 0, sum1 = 0;
			if(writedata<=20){
			  fprintf(data,"%s",filename);
			}
			

			fprintf(result,"%s",filename);
			for (i = 0; i < n; i++) {       
			  x[i]=main0(argc-1, argvs,i,data,result,writedata);
			}

			for (i = 0; i < n; i++) {
			  sum = sum + x[i];
			}
			average = sum / (float)n;

			for (i = 0; i < n; i++) {
			  sum1 = sum1 + pow((x[i] - average), 2);
			}
			variance = sum1 / (float)n;
			std_deviation = sqrt(variance);
			printf("Average = %.3f\n", average);
			printf("Standard deviation = %.4f\n", std_deviation);
	        
			fprintf(result,"\t%.3f\t%.4f\n",
				average,std_deviation);
	        
		 //        if(writedata<=20){
			//   fprintf(data,"%f\t%f\n",average,std_deviation);
			// }
			time=clock()-time;
			double seconds=((double)time)/CLOCKS_PER_SEC;
			printf("takes: %f\n",seconds);
			//finish this RNA
		}
        if(feof(r)){
            break ;
        }
        fclose(result);
    }
    fclose(r);
    fclose(result);
    fclose(not0);
    fclose(data);
    totaltime=clock()-totaltime;
    double seconds=((double)totaltime)/CLOCKS_PER_SEC;
    printf("total takes: %f\n",seconds);
  return 0;
}

/*input first the fasta file, then optionally the sample_1000.out file run on the fasta, then options*/
int main0(int argc, char *argv[],int flag,FILE* data,FILE* re,int writedata) {
  int i,input = 0, gtargs = 9;
  char **args = NULL, *name;
  HASHTBL *deleteHash;
  FILE *fp;
  Set *set;
  Options *opt;

  if (argc < 2 || !strcmp(argv[1],"--help")) {
/*print out list of options
    fprintf(stderr,"Not enough arguments\n");
    exit(EXIT_FAILURE);*/
	puts("Error: Missing input file.");
	puts("Usage: RNAprofile [OPTIONS] ... FILE");
	puts("\tFILE is an RNA sequence file containing only the sequence or in FASTA format.\n");
	print_options();
	puts("\nEXAMPLES");
	puts("1. Profile 1,000 structures sampled with gtboltzmann, default options");
	puts("\t./RNAprofile <seq_file>");
	puts("2. Profile input structures obtained from gtboltzmann with verbose output");
	puts("\t./RNAprofile -e <output samples file> -v <seq_file>\n");
	exit(EXIT_FAILURE);
  } 
  set = make_Set((char*)"output.samples");
  /* set = make_Set(argv[2]); */
  opt = set->opt;
  args = (char**)malloc(sizeof(char*)*16);
  /* Set default options for gtboltzmann */
  args[1] = (char*)"--paramdir";
  args[2] = (char*)"./data/GTparams/";
  //args[2] = (char*)"../Desktop/share/gtfold/Turner99/";
  args[3] = (char*)"-o";
  args[4] = (char*)"output";
  args[5] = (char*)"-s";
  args[6] = (char*)"1000";
  args[7] = (char*)"--scale";
  args[8] = (char*)"0.0";
  for (i = 1; i < argc-1; i++) {
    //printf("argv[%d] is %s\n",i,argv[i]);
    if (!strcmp(argv[i],"-e")) {
      if (i + 1 <= argc - 2) {
	set->structfile = argv[i+1];
	i++;
	input = 1;
      }
    }
    if (!strcmp(argv[i],"-sfold")) {
      if (i + 1 <= argc - 2) {
	set->structfile = argv[i+1];
	i++;
	input = 1;
	opt->SFOLD = 1;
      }
    }
    else if (!strcmp(argv[i],"-h")) {
      if ((i + 1 <= argc - 2) && sscanf(argv[i+1],"%lf",&(opt->HC_FREQ))) {
	opt->HC_FREQ = atof(argv[i+1]);
	if (opt->HC_FREQ < 0 || opt->HC_FREQ > 100) {
	  fprintf(stderr,"Error: invalid input %f for frequency threshold\n",opt->HC_FREQ);
	  opt->HC_FREQ = -1;
	}
	i++;
      }
    }
    else if (!strcmp(argv[i],"-p")) {
      if ((i + 1 <= argc - 2) && sscanf(argv[i+1],"%lf",&(opt->PROF_FREQ))) {
	opt->PROF_FREQ = atof(argv[i+1]);
	if (opt->PROF_FREQ < 0 || opt->PROF_FREQ > 100) {
	  fprintf(stderr,"Error: invalid input %lf for frequency threshold\n",opt->PROF_FREQ);
	  opt->PROF_FREQ = -1;
	}
	i++;
      }
    }
    else if (!strcmp(argv[i],"-c")) {
      if ((i + 1 <= argc - 2) && sscanf(argv[i+1],"%lf",&(opt->COVERAGE))) {
	opt->COVERAGE = atof(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-f")) {
      if ((i + 1 <= argc - 2) && sscanf(argv[i+1],"%d",&(opt->NUM_FHC))) {
	opt->NUM_FHC = atoi(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-s")) {
      if ((i + 1 <= argc - 2) && sscanf(argv[i+1],"%d",&(opt->NUM_SPROF))) {
	opt->NUM_SPROF = atoi(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-l")) {
      if ((i + 1 <= argc - 2) && sscanf(argv[i+1],"%d",&(opt->MIN_HEL_LEN))) {
	opt->MIN_HEL_LEN = atoi(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-u")) {
      if ((i + 1 <= argc - 2) && sscanf(argv[i+1],"%d",&(opt->NUMSTRUCTS))) {
	opt->NUMSTRUCTS = atoi(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-m")) {
      if (i + 1 <= argc - 2) {
	opt->PNOISE = atoi(argv[i+1]);
	i++;
      }
    }
    else if (!strcmp(argv[i],"-o")) {
      if (i + 1 <= argc - 2) {
	opt->OUTPUT = argv[i+1];
	name = mystrdup(argv[i+1]);
	set->structfile = strcat(name,".samples");
	args[3] = argv[i];
	args[4] = argv[i+1];
	i++;
      }
    }
    else if (!strcmp(argv[i],"-i")) {
      if (i + 1 <= argc - 2) {
	opt->INPUT = argv[i+1];
	i++;
      }
    }
    else if (!strcmp(argv[i],"-n")) {
      if (i + 1 <= argc - 2) {
	opt->NATIVE = argv[i+1];
	i++;
      }
    }
/*
    else if (!strcmp(argv[i],"-k")) {
      if (i + 1 <= argc - 2) {
	opt->CYCLES = argv[i+1];
	i++;
      }
    }
*/

/*
GTBOLTZMANN OPTIONS

    else if (!strcmp(argv[i],"-d" || !strcmp(argv[i],"--dangle"))) {
      if (i + 1 <= argc - 2) {
	if (atoi(argv[i+1]) == 0 || atoi(argv[i+1]) == 2) {
	  args[gtargs] = argv[i];
	  args[gtargs+1] = atoi(argv[i+1]);
	  gtargs += 2;
	} else {
	  fprintf(STDERR,"Wrong arguments to -d/--dangle option: ignoring option\n");
	}
	i++;
      }
    }
*/
    else if (!strcmp(argv[i],"--paramdir")) {
      if (i + 1 <= argc - 2) {
	args[1] = argv[i];
	args[2] = argv[i+1];
	i++; 
      }
    }
    else if (!strcmp(argv[i],"--limitcd")) {
      if (i+1 <= argc-2) {
	args[gtargs++] = argv[i++];
	args[gtargs++] = argv[i];
      }
    }
    else if (!strcmp(argv[i],"--useSHAPE")) {
      if (i+1 <= argc-2) {
	args[gtargs++] = argv[i++];
	args[gtargs++] = argv[i];
      }
    }
    else if (!strcmp(argv[i],"--sample")) {
      if (i+1 <= argc-2) {
	args[5] = argv[i++];
	args[6] = argv[i++];
      }
    }
    else if (!strcmp(argv[i],"-w") || !strcmp(argv[i],"--workdir")) {
      if (i+1 <= argc-2) {
	args[gtargs++] = argv[i++];
	args[gtargs++] = argv[i];
      }
    }
    else if (!strcmp(argv[i],"-v"))
      opt->VERBOSE = 1;
    else if (!strcmp(argv[i],"-g"))
      opt->GRAPH = 0;
    else if (!strcmp(argv[i],"-r"))
      opt->REP_STRUCT = 1;
    else if (!strcmp(argv[i],"-t"))
      opt->TOPDOWN = 1;
    else if (!strcmp(argv[i],"-a"))
      opt->ALTTHRESH = 0;
  }
  if (!input) {
	args[0] = (char*)"gtboltzmann";
	args[gtargs] = argv[argc-1];
	/*	for(int index =0; index<gtargs+1;index++){
	  printf("%s\n",args[index]);
	  }*/	
	boltzmann_main(gtargs+1,args);
        
  }

  input_seq(set,argv[argc-1]);
  if (opt->SFOLD)
    process_structs_sfold(set);
  else
    process_structs(set);
  reorder_helices(set);
  // print_all_helices(set);
  //printf("Total number of helix classes: %d\n",set->hc_num);

  if (set->opt->NUM_FHC)
    set->opt->HC_FREQ = set_num_fhc(set);
  else if (set->opt->HC_FREQ==-1) 
    set->opt->HC_FREQ = set_threshold_entropy(set);
  
  if (set->opt->VERBOSE) {
    printf("Threshold to find frequent helices: %.1lf%%\n",set->opt->HC_FREQ);
    printf("Number of structures processed: %d\n",set->opt->NUMSTRUCTS);
  }

  //find_bools(set);
  find_freq(set);
   int result=set->num_fhc;
     // printf("%d",result);
  
    HC **helices=set->helices;
 //    if(writedata<=20){
 //      for(int i=0;i<12;i++){
	// // fprintf(data,"%d\t",helices[i]->freq);
 //      }
 //      // fprintf(data,"%d\n",helices[12]->freq);
 //    }
   
  //   if(flag==0){
  //     //HC **helices=set->helices;
  //      for(int i=0;i<result+5;i++){
	 // fprintf(re,"%d\t",helices[i]->freq);
  //      }
  //      fprintf(re,"%d\n",helices[result+5]->freq);
  //      }
    /*
  //printf("Total number of featured helix classes: %d\n",set->num_fhc);
    if (opt->SFOLD) 
    make_profiles_sfold(set);
  else
    make_profiles(set);
    //  printf("Total number of profiles: %d\n",set->prof_num);
  //print_meta(set);
  print_profiles(set);
  if (set->opt->NUM_SPROF)
    set->opt->PROF_FREQ = set_num_sprof(set);
  else if (set->opt->PROF_FREQ == -1) {
    //set->opt->PROF_FREQ = set_p_threshold(set,P_START);
    set->opt->PROF_FREQ = set_p_threshold_entropy(set);
  }
  if (set->opt->VERBOSE) {
    if (set->opt->PROF_FREQ == -2) {
      printf("No threshold set; every profile has frequency of 1\n");
    } else
      printf("setting p to %.1f\n",set->opt->PROF_FREQ);
  }
  select_profiles(set);
  //printf("Total number of selected profiles: %d\n",set->num_sprof);
  int result=set->num_sprof;
  printf("%d",result);
  Profile **profiles=set->profiles;
  int p=set->prof_num;
    if(writedata<=20){
      for(int i=0;i<p-1;i++){
	fprintf(data,"%d\t",profiles[i]->freq);
      }
      fprintf(data,"%d\n",profiles[p-1]->freq);
    }
   
    if(flag==0){
       for(int i=0;i<p-1;i++){
	 fprintf(re,"%d\t",profiles[i]->freq);
       }
       fprintf(re,"%d\n",profiles[p-1]->freq);
       }
  if (set->opt->INPUT)
    process_one_input(set);
  if (set->opt->REP_STRUCT) {
    find_consensus(set);
    print_consensus(set);
  }
  
  if (set->opt->GRAPH) {
    fp = fopen(set->opt->OUTPUT,"w");
    init_graph(fp,set);
    i = initialize(set);
    if (set->opt->INPUT)
      print_input(fp,set);
    find_LCAs(fp,set,i);
    calc_gfreq(fp,set);
    //printGraph();
    deleteHash = MemoryDFS(set->graph);
    removeEdges(deleteHash);
    //start_trans_reductn(set->graph);
    //printGraph();
    print_edges(fp,set);
    fputs("}",fp);
    fclose(fp);
    hashtbl_destroy(deleteHash);
    }*/
  return result;
}
