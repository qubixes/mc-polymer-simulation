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
