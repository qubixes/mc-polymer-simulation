#include "lowm_modes.h"
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char** argv){
	SimProperties sp;
	
	FILE* pFile, *pFileDest;
	char* dir, destDir[10000], file[10000], exec[10000], fileDest[10000];
	if(argc<2){
		printf("Need basedir\n");
		exit(0);
	}
	dir = argv[1];
	SetSimProps(&sp, dir);
	
	sprintf(destDir, "%s/new_dat", dir);
	sprintf(exec, "mkdir -p %s", destDir);
	printf("Exec: %s\n", exec);
	system(exec);

	char* str = malloc(sizeof(char)*(sp.polSize+1));
	for(int iDev=0; iDev<sp.nDev; iDev++){
		for(long t=0; t<sp.nTime; t++){
			sprintf(file, "%s/t=%li_dev=%i.res", destDir, t*sp.dT, iDev);
			pFile = fopen(file, "w");
			fprintf(pFile, "LT= %i\n", LT);
			fprintf(pFile, "LU= %i\n", LT);
			fprintf(pFile, "LV= %i\n", LT);
			fprintf(pFile, "np= %i\n", sp.nPol);
			fprintf(pFile, "maxPolLength= %li\n", sp.polSize+1);
			fclose(pFile);
		}
	}
// 	exit(0);
	for(int iDev=0; iDev<sp.nDev; iDev++){
		for(int i=0; i<sp.nPol; i++){
			sprintf(file, "%s/ptl/pol=%i_dev=%i.res", dir, i, iDev);
			pFile = fopen(file, "r");
			if(!pFile){ printf("Error openening file %s\n", file); exit(192); }
			fscanf(pFile, "%*s %*s");
			fscanf(pFile, "%*s %*s");
			fscanf(pFile, "%*s %*s");
			for(int t=0; t<sp.nTime; t++){
				sprintf(fileDest, "%s/t=%li_dev=%i.res", destDir, t*sp.dT, iDev);
				pFileDest = fopen(fileDest, "a");
				int t,u,v;
				fscanf(pFile, "%i %i %i", &t, &u, &v);
				fscanf(pFile, "%s", str);
				fprintf(pFileDest, "len= %li\n%i %i %i\n%s\n", sp.polSize, t, u, v, str);
				fclose(pFileDest);
			}
			fclose(pFile);
			printf("\b \b\b \b\b \b\b \rpolymer %i/%i", i, sp.nPol);
			fflush(NULL);
		}
	}
	printf("\n");
}