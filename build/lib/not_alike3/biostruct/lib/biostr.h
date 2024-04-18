#include "linked_list.h"

#ifndef _BIOSTR_H_
#define _BIOSTR_H_

#define MAX_MEM_ALLOC_SEQ 1048576
#define MAX_NUM_ITEMS_ALLOC 64

struct Biostr
{
	char* seq;
	char** ids;
	int* starts;
	int* ends;
	int* hide;
	int lhand;
	int rhand;
	int num_items;
	int memAllocSeq;
	int numItemsAlloc;
	int size;
};

void Biostr_init(struct Biostr* bs);
void Biostr_init_custom(struct Biostr* bs, long long unsigned int size);
void Biostr_load_file(struct Biostr* bs, const char* input_file);
struct Biostr* Biostr_slice(struct Biostr* bs, int size, int step);
void Biostr_write_to_file(struct Biostr* bs, char* output_file);
void Biostr_filter(struct Biostr* bs, struct linkedList* headers);
void Biostr_free(struct Biostr* bs);
void Biostr_print(struct Biostr* bs);
void Biostr_sample(struct Biostr* bs, int perc, char* outfile);

#endif
