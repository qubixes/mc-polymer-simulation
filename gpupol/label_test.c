#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void PrintBit(long bitArr, long nBit){
	
	for(long i=nBit-1; i>=0; i--){
		if(bitArr&(1<<i)) printf("1");
		else printf("0");
	}
}

int main(int argc, char** argv){
	
	long maxArray=32;
	long maxLabelBit=6;
	
	long* bitArray = malloc(sizeof(long)*maxArray);
	
	for(long i=0; i<maxArray; i++) bitArray[i]=0;
	
	for(long polSize=3; polSize<maxArray; polSize++){
		long bestNumLabels=0;
		long bestLabel=0;
		long bestLabelBit=0;
		for(long nLabelBit=2; nLabelBit<=maxLabelBit && nLabelBit<polSize; nLabelBit++){
		long labelMask = (1<<nLabelBit)-1;
			for(long label=0; label< (1<<nLabelBit); label++){
				long nLabels=0;
				for(long polLabel=0; polLabel<(1<<(polSize-nLabelBit)); polLabel++){
					long totArray = (label>>1) | (polLabel<<(nLabelBit-1)) | (label<<(polSize-1));
// 					prlongf("totArray=0x%x\n", totArray);
					long labelOk=1;
// 					PrlongBit(totArray, polSize+nLabelBit-2);
					for(long i=0; i<polSize-1; i++){
						if((totArray&labelMask) == label){
							labelOk=0;
							break;
						}
						totArray >>=1;
					}
// 					prlongf("  Ok=%i\n", labelOk);
					if(labelOk) nLabels++;
				}
				if(nLabels>bestNumLabels){
					bestLabel = label;
					bestNumLabels = nLabels;
					bestLabelBit = nLabelBit;
				}
			}
		}
		printf("label=0x%lx [%li bits], polSize=%li, Number of labels=%li -> max box size = %.1lf\n", bestLabel, bestLabelBit, polSize, bestNumLabels, pow(bestNumLabels*polSize, 1/3.));
	}

}
	
	