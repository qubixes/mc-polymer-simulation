#include "file.h"

int NeedsUpdate(char* src, char* dst){
	struct stat srcAttr, dstAttr;
	if(stat(src, &srcAttr)==-1){
// 		printf("%s does not exist\n", src);
		return 0;
	}
	if(stat(dst, &dstAttr)==-1){
// 		printf("%s does not exist\n", dst);
		return 1;
	}
	
	if(srcAttr.st_mtime > dstAttr.st_mtime){
// 		printf("%s newer than %s\n", src, dst);
		return 1;
	}
	else{
// 		printf("%s older than %s\n", src, dst);
		return 0;
	}
}

int NeedsUpdatePath(char* src, char* dst, char* path){
	char* srcPath = malloc(sizeof(char)*(strlen(src)+strlen(path)+2));
	char* dstPath = malloc(sizeof(char)*(strlen(dst)+strlen(path)+2));
	
	sprintf(srcPath, "%s/%s", path, src);
	sprintf(dstPath, "%s/%s", path, dst);
	int ret = NeedsUpdate(srcPath, dstPath);
	free(srcPath); free(dstPath);
	return ret;
}

long GetTFromName(char* file){
	size_t len = strlen(file);
	char* tString = (char*) malloc(len+1);
	strcpy(tString, file);
	char *p=tString;
	long t=0;
	
	p = tString+len;
	while(*p != '/' && p>tString){
		p--;
		if(*p == '_')
			*p = '\0';
	}
	while(*p != '\0' && !(*p == 't' && *(p+1)=='=') ) p++;
	sscanf(p+2, "%li", &t);
	free(tString);
	return t;
}


int GetLastT(char* dir){
	long highestT=-1;
	char exec[3000];
	char buf[1000];
	
	sprintf(exec, "ls %s | grep 't='", dir);
	FILE* pPipe = popen(exec, "r");
	if(!pPipe){
		printf("Error opening pipe to find last file\n");
		printf("-1");
	}
	
	while(!feof(pPipe)){
		if(fscanf(pPipe, "%s", buf)<=0) continue;
		long t = GetTFromName(buf);
		if(t>highestT){
			highestT = t;
		}
	}
	return highestT;
}