#define IS_MAIN
#include "pascal_contact_map.h"

int main(int argc, char** argv){
	char* infoFile, *configFile, *outFile;
	unsigned int seed = 102874012;
	long nSamples = (long)1e8;
	if(argc<4){
		printf("Need at least three arguments: genome file, conformation file, output file\n");
		exit(192);
	}
	
	infoFile = argv[1];
	configFile = argv[2];
	outFile = argv[3];
	
	if(argc > 4)
		seed = atoi(argv[4]);
	if(argc > 5)
		nSamples = (long)atof(argv[5]);
	unsigned int rng[4];
	Seed(rng, seed);
	
	PasData* pData = PascalDataInit(infoFile, configFile);
	PCMatrix* pcm  = GeneratePCMatrix(pData);
	PrintContactMatrix(pcm, pData, outFile, rng, nSamples);
	PrintPasData(pData);
}

double ComputePC(PasData* pData, int iPol, int iMono, int jPol, int jMono){
	
	double* iPos = pData->pos[iPol][iMono];
	double* jPos = pData->pos[jPol][jMono];
	
	double rsq=0;
	
	for(int k=0; k<3; k++){
		rsq += (iPos[k]-jPos[k])*(iPos[k]-jPos[k]);
	}
	
	if(rsq>1) return 0;
	
	double pc = exp(-9./2.*rsq);
	
	return pc;
}

void PrintContactMatrix(PCMatrix* pcm, PasData* pData, char* file, unsigned int rng[4], long nSamples){
	int dAlloc = 1000;
	int cAlloc = dAlloc;
	Contact* contacts = malloc(sizeof(Contact)*cAlloc);
	long nContacts=0;
	
	for(int xBin=0; xBin<pcm->nBins[0]; xBin++){
		for(int yBin=0; yBin<pcm->nBins[1]; yBin++){
			for(int zBin=0; zBin<pcm->nBins[2]; zBin++){
				int bin = xBin + yBin*pcm->nBins[0] + zBin*pcm->nBins[0]*pcm->nBins[1];
				Node* node = pcm->matrix[bin];
				while(node){
					for(int dxBin=-1; dxBin<=1; dxBin++){
						int totXBin = dxBin+xBin;
						if(totXBin <0 || totXBin >= pcm->nBins[0]) continue;
						for(int dyBin=-1; dyBin<=1; dyBin++){
							int totYBin = dyBin+yBin;
							if(totYBin <0 || totYBin >= pcm->nBins[1]) continue;
							for(int dzBin=-1; dzBin<=1; dzBin++){
								int totZBin = dzBin+zBin;
								if(totZBin <0 || totZBin >= pcm->nBins[2]) continue;
								
								int newBin = totXBin + totYBin*pcm->nBins[0] + totZBin*pcm->nBins[0]*pcm->nBins[1];
								Node* newNode = pcm->matrix[newBin];
								while(newNode){
									if(newNode->iPol > node->iPol){
										newNode = newNode->next;
										continue;
									}
									if(newNode->iPol == node->iPol && newNode->iMono >= node->iMono){
										newNode = newNode->next;
										continue;
									}
									
									double strength = ComputePC(pData, node->iPol, node->iMono, newNode->iPol, newNode->iMono);
									if(DRng(rng)<strength){
										contacts[nContacts].iPol = node->iPol;
										contacts[nContacts].iMono = node->iMono;
										contacts[nContacts].jPol = newNode->iPol;
										contacts[nContacts].jMono = newNode->iMono;
										contacts[nContacts].strength = strength;
										nContacts++;
										
										if(cAlloc <= nContacts){
											cAlloc += dAlloc;
											contacts = realloc(contacts, sizeof(Contact)*cAlloc);
										}
									}
									newNode = newNode->next;
								}
							}
						}
					}
					node = node->next;
				}
			}
		}
	}
	
	
		
	FILE* pFile = fopen(file, "w");
	if(!pFile){
		printf("Error opening file %s\n", file);
		exit(192);
	}
	
	fprintf(pFile, "#nPol= %i\n", pData->nPol);
	fprintf(pFile, "#maxLen= %i\n", pData->maxNMono);
	fprintf(pFile, "#nContacts = %li\n", nContacts);
	for(int iPol=0; iPol<pData->nPol; iPol++)
		fprintf(pFile, "lin %i\n", pData->nMono[iPol]);
	
	for(int iContact=0; iContact<nContacts && iContact < nSamples; iContact++){
		int j = DRng(rng)*(nContacts-iContact);
		fprintf(pFile, "%i %i %i %i %lf\n", contacts[j].iPol, contacts[j].iMono, contacts[j].jPol, contacts[j].jMono, 1.0);
		contacts[j] = contacts[nContacts-iContact-1];
	}
	fclose(pFile);
}

int Index(PCMatrix* pcm, double xyz[3]){
	int index=0;
	
	int mult=1;
	for(int k=0; k<3; k++){
		int xyzIndex = (int)((xyz[k]-pcm->xyzStart[k])/pcm->dxyz);
		index += mult*xyzIndex;
		mult *= pcm->nBins[k];
	}
	return index;
}

Node* NewNode(int iPol, int iMono){
	Node* node = malloc(sizeof(Node));
	node->iPol = iPol;
	node->iMono = iMono;
	return node;
}

PCMatrix* GeneratePCMatrix(PasData* pData){
	PCMatrix* pcm = malloc(sizeof(PCMatrix));
	
// 	pcm->nBinXYZ = MAX(1, (int)(0.5*pow(pData->nTotMono, 1./3.)+0.5));
	
	/// Choose the cut-off as the box sizes. Then you only need to search the box+8 adjacent ones.
	pcm->dxyz = 1;
	
	
	double minXYZ[3] = {1e9,1e9,1e9}, maxXYZ[3]={-1e9,-1e9,-1e9};
	for(int iPol=0; iPol< pData->nPol; iPol++){
		for(int iMono=0; iMono<pData->nMono[iPol]; iMono++){
			for(int k=0; k<3; k++){
				maxXYZ[k] = MAX(pData->pos[iPol][iMono][k], maxXYZ[k]);
				minXYZ[k] = MIN(pData->pos[iPol][iMono][k], minXYZ[k]);
			}
		}
	}
	
	pcm->nTotBins=1;
	for(int k=0; k<3; k++){
		pcm->nBins[k] = (maxXYZ[k]-minXYZ[k])/pcm->dxyz + 1;
		pcm->nTotBins *= pcm->nBins[k];
		pcm->xyzStart[k] = minXYZ[k];
	}
	
	pcm->matrix = malloc(sizeof(Node*)*pcm->nTotBins);
	for(int i=0; i<pcm->nTotBins; i++) pcm->matrix[i] = NULL;
	
	for(int iPol=0; iPol<pData->nPol; iPol++){
		for(int iMono=0; iMono<pData->nMono[iPol]; iMono++){
			int index = Index(pcm, pData->pos[iPol][iMono]);
			Node* newNode = NewNode(iPol, iMono);
			newNode->next = pcm->matrix[index];
			pcm->matrix[index] = newNode;
		}
	}
	return pcm;
}

