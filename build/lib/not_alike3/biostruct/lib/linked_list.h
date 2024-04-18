#ifndef _LINKLIST_H_
#define _LINKLIST_H_

struct Node
{
	char* header;
	struct Node* next;
};

struct linkedList
{
	struct Node* head;
	struct Node* tail;
};

struct linkedList* create_linked_list(void);
struct Node* create_node(const char header[]);
void add_node(struct linkedList* ll, const char header[]);
void free_linked_list(struct linkedList* ll);
void free_node(struct Node* n);
void insert_at_position(struct linkedList* ll, const char header[], int position);
void print_linked_list(struct linkedList* ll);
void load_to_linkedList(struct linkedList* ll, const char filename[]);

#endif
