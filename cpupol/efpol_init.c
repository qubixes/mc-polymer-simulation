#include "efpol.h"

int GetDtMax(char* dir, long* dt, long* tMax){
	char exec[10000];
	char* file = malloc(10000*sizeof(char));
// 	char* retFile = malloc(10000*sizeof(char));
	FILE* pPipe;
	int i;
	int fileFound=0;
	long t;
	long minT=(long)1e13, maxT=(long)0, secMinT=(long)1e13;
	
	sprintf(exec, "ls %s | grep 't='", dir);
	pPipe = popen(exec, "r");
// 	printf("Command: %s\n", exec);
	if(!pPipe){
		printf("error opening pipe\n");
		exit(0);
	}
	while(fscanf(pPipe, "%s", file)>0){
		i=0; 
		while(file[i] != '_' && file[i] != '\0') i++;
		if(file[i] != '_') continue;
		file[i]='\0';
		t = atol(file+2);
		file[i]='_';
// 		printf("t = %li\n", t);
		if(t<minT){
			secMinT = minT;
			minT = t;
		}
		else if(t<secMinT)
			secMinT = t;
		
		if(t>maxT)
			maxT = t;
		fileFound=1;
	}
	pclose(pPipe);
	
	if(!fileFound) return 0;
	
	*dt = secMinT-minT;
	*tMax = maxT;
	
	if(secMinT == (long)1e13 || maxT ==0){
		*dt=0;
		*tMax=maxT;
	}
// 	printf("tMax = %li\n", maxT);
// 	exit(0);
	return 1;
}

void CSInit(CurState* cs, int BL, int polSize, int nPol, char* dir){
	cs->polSize = polSize;
	cs->nPol = nPol;
	cs->intsPerPol = cs->polSize/8+((cs->polSize%8==0)?0:1);
	cs->intsPerMono = cs->nPol/8+((cs->nPol%8==0)?0:1);
	SetConstants(&cs->con, BL);

	cs->polSize = polSize;
	cs->lattice    = malloc(sizeof(uint)*cs->con.LAT_ALLOC);
	cs->coorPol    = malloc(sizeof(uint)*cs->polSize*cs->nPol);
	cs->unitPol    = malloc(sizeof(uint)*cs->intsPerMono*cs->polSize);
	cs->unit1Cache = malloc(sizeof(uint)*cs->intsPerMono);
	cs->unit2Cache = malloc(sizeof(uint)*cs->intsPerMono);
	cs->coorCache  = malloc(sizeof(uint)*cs->nPol);
	cs->allocated  = TRUE;
	if(dir){
		cs->dir        = malloc(sizeof(char)*(strlen(dir)+200));
		sprintf(cs->dir, "%s/N_%i", dir, cs->polSize);
		char exec[1000];
		sprintf(exec, "mkdir -p %s\n", cs->dir);
		system(exec);
	}
	CSClear(cs);
}

void CSClear(CurState* cs){
	for(int coor=0; coor<cs->con.LAT_ALLOC; coor++){
		cs->lattice[coor]=0;
	}
	
	for(int i=0; i<cs->polSize*cs->nPol; i++) {
		cs->coorPol[i] = 0;
	}
	for(int i=0; i<cs->intsPerMono*cs->polSize; i++)
		cs->unitPol[i] = 0;
}

void PCInit(CurState* cs, PolyConfig* pcfg){
	pcfg->xyz = malloc(sizeof(CoorI)*cs->polSize);
	pcfg->tuv = malloc(sizeof(CoorI)*cs->polSize);
}

void PrintLattice(CurState* cs){
	
	for(int t=0; t<cs->con.L; t++){
		for(int u=0; u<cs->con.L; u++){
			for(int v=0; v<cs->con.L; v++){
				int coor = TUVtoCoor(t,u,v, &cs->con);
				if(OccupLattice(coor, cs))
					printf("1");
				else
					printf("0");
			}
			printf("\n");
		}
		printf("\n");
	}
}

void DoubleLattice(CurState* cs, CurState* cs_next){
	int t,u,v;
	CSClear(cs_next);
	for(int iPol=0; iPol<cs->nPol; iPol++){
		int nUnit=0;
		for(int iMono=0; iMono<cs->polSize; iMono++){
			int coor = cs->coorPol[iPol+iMono*cs->nPol];
			t=TCoor(coor,&cs->con); u=UCoor(coor,&cs->con); v=VCoor(coor,&cs->con);
			int newCoor = TUVtoCoor(2*t,2*u,2*v,&cs_next->con);
			uint unit = (cs->unitPol[iPol/8+iMono*cs->intsPerMono]>>(4*(iPol%8)))&0xf;
			nUnit+=unit?1:0;
			cs_next->unitPol[iPol/8+(8*iMono+3)*cs_next->intsPerMono] |= unit<<(4*(iPol%8));
			cs_next->unitPol[iPol/8+(8*iMono+7)*cs_next->intsPerMono] |= unit<<(4*(iPol%8));
			SetLattice(newCoor, cs_next);
			for(int i=0; i<4; i++)
				cs_next->coorPol[iPol+(8*iMono+i)*cs->nPol]=newCoor;
			newCoor = AddUnitToCoor(unit, newCoor, &cs_next->con);
			for(int i=4; i<8; i++)
				cs_next->coorPol[iPol+(8*iMono+i)*cs->nPol]=newCoor;
			SetLattice(newCoor, cs_next);
// 			if(nUnit==2)
// 				printf("Uh oh!\n");
		}
// 		printf("nUnit=%i\n", nUnit);
		
		if(nUnit==2){
// 			CSClear(cs_next);
			int coorSave[5];
			uint unit=0x0;
			int coor=0;
			for(int iMono=0; iMono<cs->polSize; iMono++){
				unit = (cs->unitPol[iPol/8+iMono*cs->intsPerMono]>>(4*(iPol%8)))&0xf;
				if(unit && unit<8){
					coor=cs->coorPol[iPol+iMono*cs->nPol];
					break;
				}
			}
			int try=0;
			t=TCoor(coor,&cs->con);u=UCoor(coor,&cs->con); v=VCoor(coor,&cs->con);
			int baseCoor = TUVtoCoor(2*t,2*u,2*v,&cs_next->con);
			int combo = (~unit)&0xf;
			int fail=1;
			for(try=0; try<4 && fail; try++){
// 				if(try>0)
// 					printf("Attempt no. %i (pol=%i)\n", try, iPol);
				
				///First find the unit vector.
				uint dst1= tab.transTable[combo*NOPT+try]&0xf;
				uint dst2= combo;
				uint dst3= (tab.transTable[combo*NOPT+try]>>4)&0xf;
				
				uint newUnits[5]={unit,unit, dst1, dst2, dst3};
				
				int iMono=0;
	// 			printf("newCoor=%x\n", newCoor);
				int newCoor=baseCoor;
				for(int i=0; i<5; i++){
	// 				printf("unit=%x\n", newUnits[i]);
					for(int j=0; j<cs_next->polSize/5+((i<cs_next->polSize%5)?1:0)-1; j++){
						cs_next->coorPol[iPol+iMono*cs_next->nPol] = newCoor;
						cs_next->unitPol[iPol/8+iMono*cs_next->intsPerMono] &= ~(0xf<<(4*(iPol%8)));
						iMono++;
					}
					cs_next->unitPol[iPol/8+iMono*cs_next->intsPerMono] &= ~(0xf<<(4*(iPol%8)));
					cs_next->unitPol[iPol/8+iMono*cs_next->intsPerMono] |= newUnits[i]<<(4*(iPol%8));
					cs_next->coorPol[iPol+iMono*cs_next->nPol] = newCoor;
					iMono++;
					newCoor = AddUnitToCoor(newUnits[i], newCoor, &cs_next->con);
					if(i>=2 && i<4 && OccupLattice(newCoor, cs_next)){
// 						printf("\nRemoving: unit=%i, coor=%x, dst1=%x, dst2=%x, dst3=%x, i=%i, iPol=%i\n", unit, newCoor, dst1, dst2, dst3,i, iPol);

						for(int j=2; j<i; j++){
// 							printf("Unsetting coor: %i (pId=%i)\n", coorSave[j], iPol);
							UnsetLattice(coorSave[j], cs_next);
// 						}
						}
						break;
					}
// 						printf("\nUhoh: newCoor=%x!\n", newCoor);
// 						exit(0);
					SetLattice(newCoor, cs_next);
					coorSave[i]=newCoor;
					if(i==4) fail=0;
				}
			}
			if(fail){
				printf("Failed to find an empty spot...\n");
				exit(0);
			}
		}
	}
// 	PrintPol(95,cs);
// 	PrintLattice(cs);
}

void SetConstants(Constants* con, int BL){
	con->BL          = BL;
	con->L           = 1<<BL;
	con->ALLOC_T     = 2*con->L;
	con->ALLOC_U     = 2*con->L;
	con->ALLOC_V     = con->L;
	con->LAT_SIZE    = (con->L*con->L*con->L);
	con->LAT_ALLOC   = (con->ALLOC_T*con->ALLOC_U*con->ALLOC_V)/32;
	con->LARGE_T_UNIT= 1;
	con->LARGE_U_UNIT= con->ALLOC_T;
	con->LARGE_V_UNIT= con->ALLOC_T*con->ALLOC_U;
	con->LARGE_W_UNIT= con->ALLOC_T*con->ALLOC_U*con->ALLOC_V*2;
	con->LARGE_TUVW  = con->LARGE_T_UNIT|con->LARGE_U_UNIT|con->LARGE_V_UNIT|con->LARGE_W_UNIT;
	con->ADD_BASE    = con->L+con->L*con->ALLOC_T+con->L*con->ALLOC_T*con->ALLOC_U;
	con->T_MASK      = (con->ALLOC_T-1);
	con->U_MASK      = (con->ALLOC_U-1)*con->ALLOC_T;
	con->V_MASK      = (con->ALLOC_V-1)*con->ALLOC_T*con->ALLOC_U;
	con->BOUND_MASK  = (con->L-1)|((con->L-1)*con->ALLOC_T)|((con->L-1)*con->ALLOC_T*con->ALLOC_U);
}

void ResultInit(Config* cfg, CurState* cs, Result* res){
	PCInit(cs, &res->pcfg);
	res->gen = malloc(sizeof(double)*cs->polSize);
	res->genCounts = malloc(sizeof(int)*cs->polSize);
	for(int i=0; i<cs->polSize; i++){
		res->gen[i]=0;
		res->genCounts[i]=0;
	}
}

void CreatePolymers(Config* cfg, CurState* cs){
	int startCoor, endCoor;
	int curT, curU, curV, delta, nTot;
	int nTPol, nUPol, nVPol;
	
	int LT=cs->con.L, LU=cs->con.L, LV=cs->con.L;
	
	for(delta=1; delta<LT && delta<LU && delta<LV; delta++){
		nTot = 1;
		nTot *= LT/delta;
		nTot *= LU/delta;
		nTot *= LV/(delta+1);
		if(nTot < cs->nPol){
			delta--;
			break;
		}
	}
// 	printf("delta = %i\n", delta);
	nTPol = LT/delta;
	nUPol = LU/delta;
	nVPol = LV/(delta+1);
	for(int iPol=0; iPol<cs->nPol; iPol++){
		curT = delta*(iPol%nTPol);
		curV = (delta+1)*((iPol/nTPol)%nVPol);
		curU = delta*(iPol/nTPol/nVPol);
		startCoor = TUVtoCoor(curT, curU, curV,   &cs->con);
		endCoor   = TUVtoCoor(curT, curU, curV+1, &cs->con);
		
		cs->coorPol[iPol] = startCoor;
// 		PrintCoor(startCoor, &cs->con);printf("-");
// 		PrintCoor(endCoor, &cs->con);printf("\n");
		cs->unitPol[iPol/8] |= 0x4<<((iPol%8)*4); // v-direction
		//Every other element of unitPol is already initialized to stored length (so we don't need to change it).
		SetLattice(startCoor, cs); SetLattice(endCoor, cs);
		
		for(int iMono=1; iMono<cs->polSize; iMono++){
			cs->coorPol[iPol+cs->nPol*iMono] = endCoor;
		}
#if POL_TYPE == POL_LIN
		cs->unitPol[iPol/8+cs->intsPerMono*(cs->polSize-1)] |= 0xf<<(4*(iPol%8)); //End of the polymer
#else
		cs->unitPol[iPol/8+cs->intsPerMono*(cs->polSize-1)] |= 0xb<<((iPol%8)*4); //(-v) direction
#endif
	}
}

void AddTrans(uint src1, uint src2, uint dst1, uint dst2){
	int occupId = 16*NOPT*src1+NOPT*src2+tab.curNTrans[16*src1+src2];
	tab.transTable[occupId] = dst1|(dst2<<4);
	tab.accRat[occupId] = 1;
	if(dst1 || dst2){
		uint res;
		if(dst1 != 0xf){
			if(!ValidateAddUnitVectors(dst1, (~src1)&0xf, &res)){
				printf("Error adding trans: [%x,%x] => [%x,%x]\n", src1, src2, dst1, dst2);
				exit(0);
			}
		}
		else{
			if(!ValidateAddUnitVectors(src2, (~dst2)&0xf, &res)){
				printf("Error adding trans: [%x,%x] => [%x,%x]\n", src1, src2, dst1, dst2);
				exit(0);
			}
		}
// 		if(XUnit(res) == -1)
// 			tab.accRat[occupId] = exp(-cfg.EField);
		tab.transTable[occupId] |= res<<8;
	}
	if(src1 && src2) tab.remOccup[occupId]=1;
	else tab.remOccup[occupId] = 0;
	if(dst1 && dst2) tab.addOccup[occupId]=1;
	else tab.addOccup[occupId] = 0;
	tab.curNTrans[16*src1+src2]++;
}

void SetTransTable(Constants* con, int polModel){
	unsigned int unit1, unit2, unitRes;
	unsigned int valid;
	
	for(int i=0; i<16*16; i++) tab.curNTrans[i]=0;
	
	for(unit1=0x1; unit1<0xf; unit1++){
		if(!IsValid(unit1)) continue;
		for(unit2=unit1+1; unit2<0xf; unit2++){
			if(!IsValid(unit2)) continue;
			
			valid = ValidateAddUnitVectors(unit1, unit2, &unitRes);
			if(!valid) continue;
			AddTrans(unit1, unit2, unitRes, 0);
			AddTrans(unit1, unit2, 0, unitRes);
			AddTrans(unit2, unit1, unitRes, 0);
			AddTrans(unit2, unit1, 0, unitRes);
			
			AddTrans(unitRes, 0, unit1, unit2);
			AddTrans(unitRes, 0, unit2, unit1);
			AddTrans(0, unitRes, unit1, unit2);
			AddTrans(0, unitRes, unit2, unit1);
		}
	}
// 	unit1 = 0;
	for(unsigned unit2 = 1; unit2<0xf; unit2++){
		if(!IsValid(unit2)) continue;
		for(int i=0; i<4; i++){
			AddTrans(0, unit2, unit2, 0);
			AddTrans(unit2, 0, 0, unit2);
		}
	}
	
	uint res;
	for(uint unit1=1; unit1<0xf; unit1++){
		if(!IsValid(unit1)) continue;
		for(uint unit2=1; unit2<0xf; unit2++){
			if(!IsValid(unit2)) continue;
			if(!ValidateAddUnitVectors(unit1, (~unit2)&0xf, &res)) continue;
			for(uint unit3=1; unit3<0xf; unit3++){
				if(!IsValid(unit3)) continue;
				if(polModel == SL_QUAD){
					if(ValidateAddUnitVectors(unit1, unit3, &res)) continue;
				}
				for(uint unit4=1; unit4<0xf; unit4++){
					if(!IsValid(unit4)) continue;
					uint addUnit = SmallToLargeUnit(unit1, con) + SmallToLargeUnit(unit3, con) + SmallToLargeUnit((~unit2)&0xf, con) + SmallToLargeUnit((~unit4)&0xf, con);
					if(addUnit%con->LARGE_TUVW == 0){
						AddTrans(unit1, unit3, unit2, unit4);
						AddTrans(unit1, unit3, unit2, unit4);
						if(polModel == SL_EQUAL){
							if(ValidateAddUnitVectors(unit1, unit3, &res)){
								for(int i=0; i<2; i++){
									AddTrans(unit1, unit3, unit2, unit4);
									AddTrans(unit1, unit3, unit2, unit4);
								}
							}
						}
					}
				}
			}
		}
	}
	
	for(uint option=0; option<256; option++){
		if(tab.curNTrans[option] == 2){
			for(int i=0; i<2; i++){
				uint dst1 = tab.transTable[option*NOPT+i]&0xf;
				uint dst2 = (tab.transTable[option*NOPT+i]>>4)&0xf;
				AddTrans(option/16, option%16, dst1, dst2);
			}
		}
		if(tab.curNTrans[option] == 4){
			for(int i=0; i<4; i++){
				uint dst1 = tab.transTable[option*NOPT+i]&0xf;
				uint dst2 = (tab.transTable[option*NOPT+i]>>4)&0xf;
				AddTrans(option/16, option%16, dst1, dst2);
			}
		}
		if(tab.curNTrans[option] == 8){
			for(int i=0; i<8; i++){
				uint dst1 = tab.transTable[option*NOPT+i]&0xf;
				uint dst2 = (tab.transTable[option*NOPT+i]>>4)&0xf;
				AddTrans(option/16, option%16, dst1, dst2);
			}
		}
	}
	
	for(uint unit1=0; unit1<0xf; unit1++){
		if(!IsValid(unit1)) continue;
// 		int nMove=-1;
// 		if(polModel == SL_DOUBLE) nMove=4;
// 		else if(polModel == SL_EQUAL) nMove=2;
// 		else if(polModel == SL_HALF) nMove=1;
// 		else {printf("polModel????\n"); exit(0);}
		for(int i=0; i<polModel; i++){
			AddTrans(unit1, 0xf, 0x0, 0xf);
			AddTrans(0xf, unit1, 0xf, 0x0);
		}
		
		AddTrans(0x0, 0xf, unit1, 0xf);
		AddTrans(0xf, 0x0, 0xf, unit1);
	}
	for(uint option=0; option<256; option++){
		while(tab.curNTrans[option] < 16){
			AddTrans(option/16, option%16, 0, 0);
		}
	}
// 	PrintTrans(); 
// 	exit(0);
}

void ConfigInit(){
	char exec[1000];

	sprintf(exec, "mkdir -p %s", cfg.dir);
	system(exec);
	
	Seed(cfg.rngState, cfg.seed);
}
