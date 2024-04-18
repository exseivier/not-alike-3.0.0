#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linked_list.h"
#include "biostr.h"

  /*****************************************/
 /***	Static functions	Begin	***/
/*****************************************/

static char*
__open_file(const char* input_file)
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

static void
__resize_MemAlloc_Seq(struct Biostr* bs, int expected_size)
{
	bs->memAllocSeq = (bs->memAllocSeq + expected_size) * 2;
	bs->seq = (char*) realloc (bs->seq, sizeof(char) * bs->memAllocSeq + 1);
	bs->seq[bs->memAllocSeq] = '\0';
	return;
}

static void
__resize_NumItems_Alloc(struct Biostr* bs)
{
	bs->numItemsAlloc = bs->numItemsAlloc * 2;
	bs->ids = (char**) realloc (bs->ids, sizeof(char*) * bs->numItemsAlloc + 1);
	bs->starts = (int*) realloc (bs->starts, sizeof(int) * bs->numItemsAlloc + 1);
	bs->ends = (int*) realloc (bs->ends, sizeof(int) * bs->numItemsAlloc + 1);
	bs->hide = (int*) realloc (bs->hide, sizeof(int) * bs->numItemsAlloc + 1);
	bs->ids[bs->numItemsAlloc] = '\0';
	bs->starts[bs->numItemsAlloc] = '\0';
	bs->ends[bs->numItemsAlloc] = '\0';
	bs->hide[bs->numItemsAlloc] = '\0';
	return;
}

/*
 * Tested ok
 */
static char* __extract_sequence(struct Biostr* bs, const char header[])
{
	int i = 0;
	while (NULL != bs->ids[i])
	{
		if (strcmp(bs->ids[i], header) == 0)
		{
			char* sequence = (char*) calloc (bs->ends[i] - bs->starts[i] + 2, sizeof(char));
			for (int j = bs->starts[i]; j <= bs->ends[i]; j++)
			{
				sequence[j - bs->starts[i]] = bs->seq[j];
			}
			sequence[bs->ends[i] - bs->starts[i] + 1] = '\0';
			return sequence;
		}
		i++;
	}
	return NULL;
}

  /*****************************************/
 /***	Static functions	End	***/
/*****************************************/

void
Biostr_init(struct Biostr* bs)
{
	bs->seq = (char*) calloc (MAX_MEM_ALLOC_SEQ + 1, sizeof(char));
	bs->ids = (char**) calloc (MAX_NUM_ITEMS_ALLOC + 1, sizeof(char*));
	bs->starts = (int*) calloc (MAX_NUM_ITEMS_ALLOC + 1, sizeof(int));
	bs->ends = (int*) calloc (MAX_NUM_ITEMS_ALLOC + 1, sizeof(int));
	bs->hide = (int*) calloc (MAX_NUM_ITEMS_ALLOC + 1, sizeof(int));
	bs->lhand = 0;
	bs->rhand = 0;
	bs->size = 0;
	bs->num_items = 0;
	bs->memAllocSeq = MAX_MEM_ALLOC_SEQ;
	bs->numItemsAlloc = MAX_NUM_ITEMS_ALLOC;

	bs->seq[MAX_MEM_ALLOC_SEQ] = '\0';
	bs->ids[MAX_NUM_ITEMS_ALLOC] = '\0';
	bs->starts[MAX_NUM_ITEMS_ALLOC] = '\0';
	bs->ends[MAX_NUM_ITEMS_ALLOC] = '\0';
	bs->hide[MAX_NUM_ITEMS_ALLOC] = '\0';
	return;
}

void
Biostr_init_custom(struct Biostr* bs, unsigned long long size)
{
	bs->seq = (char*) calloc (size+1, sizeof(char));
	bs->ids = (char**) calloc (MAX_NUM_ITEMS_ALLOC + 1, sizeof(char*));
	bs->starts = (int*) calloc (MAX_NUM_ITEMS_ALLOC + 1, sizeof(int));
	bs->ends = (int*) calloc (MAX_NUM_ITEMS_ALLOC + 1, sizeof(int));
	bs->hide = (int*) calloc (MAX_NUM_ITEMS_ALLOC + 1, sizeof(int));
	bs->lhand = 0;
	bs->rhand = 0;
	bs->size = 0;
	bs->num_items = 0;
	bs->memAllocSeq = size;
	bs->numItemsAlloc = MAX_NUM_ITEMS_ALLOC;

	bs->seq[size] = '\0';
	bs->ids[MAX_NUM_ITEMS_ALLOC] = '\0';
	bs->starts[MAX_NUM_ITEMS_ALLOC] = '\0';
	bs->ends[MAX_NUM_ITEMS_ALLOC] = '\0';
	bs->hide[MAX_NUM_ITEMS_ALLOC] = '\0';
	return;
}

void
Biostr_load_file(struct Biostr* bs, const char* input_file)
{
	char* file_str = __open_file(input_file);
	char* line = strtok(file_str, "\n");
	int line_size = strlen(line);
	while (NULL != line)
	{

		if (line[0] == '>')
		{
			bs->num_items++;

			if (bs->num_items >= bs->numItemsAlloc)
				__resize_NumItems_Alloc(bs);

			bs->ids[bs->num_items-1] = (char*) calloc (line_size + 1, sizeof(char));
			for(int i = 0; i < line_size; i++)
			{
				bs->ids[bs->num_items-1][i] = line[i];
			}
			bs->ids[bs->num_items-1][line_size] = '\0';
			bs->starts[bs->num_items-1] = bs->rhand;
			bs->hide[bs->num_items-1] = 0;
		}
		else
		{
			bs->rhand = bs->lhand + line_size;

			if (bs->rhand >= bs->memAllocSeq)
				__resize_MemAlloc_Seq(bs, bs->rhand);

			for(int i = 0; i < bs->rhand; i++)
			{
				bs->seq[bs->lhand + i] = line[i];
			}
			bs->seq[bs->rhand] = '\0';
			bs->lhand = bs->rhand;
			bs->ends[bs->num_items-1] = bs->rhand-1;
			bs->size = bs->rhand;
		}
		line = strtok(NULL, "\n");
		if (NULL != line)
			line_size = strlen(line);
	}
	free(file_str);
	free(line);
	return;
}

struct Biostr*
Biostr_slice(struct Biostr* bs, int size, int step)
{
	struct Biostr* ret_bs = (struct Biostr*) calloc (1, sizeof(struct Biostr));
	Biostr_init_custom(ret_bs, bs->size);
	unsigned int start = 0;
	unsigned int end = 0;
	unsigned int frag_num = 0;
	char id[150];
	int len_id = 0;
	ret_bs->size = bs->size;
	for (int i = 0; i < bs->size; i++)
	{
		ret_bs->seq[i] = bs->seq[i];
	}
	for (int i = 0; i < bs->num_items; i++)
	{
		start = bs->starts[i];
		end = bs->ends[i];
		while (start <= end)
		{
			sprintf(id, "%s_frag_%d", bs->ids[i], frag_num);
			len_id = strlen(id);
			ret_bs->ids[frag_num] = (char*) calloc (len_id + 1, sizeof(char));
			ret_bs->ids[frag_num][len_id] = '\0';
			for (int i = 0; i < len_id; i++)
			{
				ret_bs->ids[frag_num][i] = id[i];
			}
			ret_bs->starts[frag_num] = start;
			ret_bs->ends[frag_num] = start + size - 1 < end ? start + size - 1 : end;
			ret_bs->hide[frag_num] = 0;
			start = start + size;
			frag_num++;
			ret_bs->num_items++;
		}
	}
	return ret_bs;
}

/*
 * tested ok
 */

void
Biostr_write_to_file(struct Biostr* bs, char* filename)
{
	FILE* FH = fopen(filename, "w+");
	if (NULL == FH)
	{
		printf("FileIOError: Unable to open file %s\n", filename);
		return;
	}
	int i = 0;
	char* sequence = NULL;
	while (NULL != bs->ids[i])
	{
		if (bs->hide[i] == 0)
		{
			sequence = __extract_sequence(bs, bs->ids[i]);
			fprintf(FH, "%s\n%s\n", bs->ids[i], sequence);
			free(sequence);
		}
		i++;		
	}
	fclose(FH);
	return;
}

/*
 * Tested ok
 */

void
Biostr_filter(struct Biostr* bs, struct linkedList* headers)
{
	//	TODO: Test this function.
	//	XXX: Write load_headers function.
	//	It must take a text file of BLAST hits and upload
	//	all headers to linkedList data structure.
	int i = 0;
	while (NULL != bs->ids[i])
	{
		struct Node* current = headers->head;
		struct Node* trash = NULL;
		struct Node* before = NULL;
		while (NULL != current)
		{
			if (strcmp(current->header, bs->ids[i]) == 0)
			{
				bs->hide[i] = 1;
				if (NULL != before)
				{
					trash = current;
					before->next = current->next;
					current = current->next;
					free_node(trash);
					break;
				}
				else
				{
					trash = current;
					current = current->next;
					free_node(trash);
					headers->head = current;
					break;
				}
			}
			else
			{
				before = current;
				current = current->next;
			}
		}
		i++;
	}
	return;
}

void
Biostr_free(struct Biostr* bs)
{
	int i = 0;
	while (NULL != bs->ids[i])
	{
		free(bs->ids[i]);
		i++;
	}
	free(bs->seq);
	free(bs->starts);
	free(bs->ends);
	free(bs->hide);
	free(bs);
	return;
}

void
Biostr_print(struct Biostr* bs)
{
	printf("Sequence: %s\n", bs->seq);
	printf("Size: %d\n", bs->size);
	for (int i = 0; i < bs->num_items; i++)
	{
		printf("ID: %s\tStart: %d\tEnd: %d\tHide: %d\n", \
				bs->ids[i], \
				bs->starts[i], \
				bs->ends[i], \
				bs->hide[i]);
	}
	return;
}

void
Biostr_sample(struct Biostr* bs, int perc, char* outfile)
{
	/*
	 * Samples sequences from DNA structure and writes them to file.
	 * 	perc means the percent of sequences to sample.
	 */

	int randint = 0;
	int i = 0;
	FILE* out = fopen(outfile, "w");
	srand(time(0));
	while (bs->ids[i] != NULL)
	{
		randint = rand() % 100;
		if (randint < perc)
		{
			
			char* seqCarrier = (char*)calloc(sizeof(char), (bs->end[i] - bs->start[i] + 2));
			int j = 0;
			int k = bs->start[i];
			while (k < bs->end[i]+1)
			{
				seqCarrier[j] = bs->seq[k];
				j++;
				k++;
			}
			seqCarrier[j] = '\0';
			fprintf(out, "%s\n%s\n", bs->ids[i], seqCarrier);
			free(seqCarrier);
			seqCarrier = NULL;
		}
		i++;
	}
	fclose(out);

}
