#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linked_list.h"
#include "fileio.h"

struct Node*
create_node(const char header[])
{
	struct Node* node = (struct Node*) calloc (1, sizeof(struct Node));
	int len_header = strlen(header);
	node->header = (char*) calloc (len_header + 1, sizeof(char));
	node->next = NULL;
	for (int i = 0; i < len_header; i++)
	{
		node->header[i] = header[i];
	}
	node->header[len_header] = '\0';

	return node;
}

struct linkedList*
create_linked_list(void)
{
	struct linkedList* ll = (struct linkedList*) calloc (1, sizeof(struct linkedList));
	ll->head = NULL;
	ll->tail = NULL;
	return ll;
}

void
add_node(struct linkedList* ll, const char header[])
{
	struct Node* node = create_node(header);
	if (NULL == ll->head)
	{
		ll->head = node;
		ll->tail = node;
	}
	else
	{
		ll->tail->next = node;
		ll->tail = node;
	}
}

void
free_node(struct Node* n)
{
	free(n->header);
	free(n);
	return;
}

void
free_linked_list(struct linkedList* ll)
{
	struct Node* tmp_node = ll->head;
	while (NULL != tmp_node)
	{
		free(tmp_node->header);
		tmp_node = tmp_node->next;
	}
	free(ll->head);
	free(ll->tail);
	free(ll);
}

void
insert_at_position(struct linkedList* ll, const char header[], int position)
{
	int counter = 0;
	struct Node* current = ll->head;
	struct Node* before = NULL;
	struct Node* node = create_node(header);
	while (NULL != current)
	{
		if (counter == position)
		{
			if (NULL == before)
			{
				node->next = current;
				current = node;
				return;
			}
			else
			{
				before->next = node;
				node->next = current;
				return;
			}
		}
		before = current;
		current = current->next;
		counter++;
	}
	return;
}

void
print_linked_list(struct linkedList* ll)
{
	struct Node* tmp_node = ll->head;
	int node_counter = 0;
	while (NULL != tmp_node)
	{
		printf("node %d: %s\n", node_counter++, tmp_node->header);
		tmp_node = tmp_node->next;
	}
	return;
}

void
load_to_linkedList(struct linkedList* ll, const char filename[])
{
	// TODO: write load to linkedList from file function.
	// XXX: Also test Biostr_filter, Biostr_wrtie_to_file
	// XXX: and the other functions.
	
	char* file_str = open_file(filename);
	char* line = strtok(file_str, "\n");
	while (NULL != line)
	{
		add_node(ll, line);
		line = strtok(NULL, "\n");
	}
	return;
}
