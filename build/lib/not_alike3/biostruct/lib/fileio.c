
#include <stdio.h>
#include <stdlib.h>

char*
open_file(const char* input_file)
{
	FILE* FH = fopen(input_file, "rb");
	fseek(FH, 0L, SEEK_END);
	long long unsigned int file_size = ftell(FH);
	rewind(FH);
	char* file_str = (char*) calloc (file_size + 1, sizeof(char));
	if (NULL == file_str)
		return NULL;
	fread(file_str, sizeof(char), file_size, FH);
	fclose(FH);
	return file_str;
}


