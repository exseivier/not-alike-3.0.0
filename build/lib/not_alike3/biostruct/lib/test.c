#include <stdio.h>
#include <stdlib.h>
#include "biostr.h"
#include "linked_list.h"

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		printf("ArgError: No arguments!\n");
		return 1;
	}
	struct Biostr* seqs = (struct Biostr*) calloc (1, sizeof(struct Biostr));
	Biostr_init(seqs);
	Biostr_load_file(seqs, argv[1]);
	Biostr_print(seqs);
	struct Biostr* sliced_seqs = Biostr_slice(seqs, 6, 2);
	printf("\nSliced\n");
	Biostr_print(sliced_seqs);
	Biostr_free(seqs);
	printf("-----------------------------------------------------\n");
	printf("Testing Biostr_wrtie_to_file and Biostr_filter\n");
	printf("Creating linkedList ... OK\n");
	struct linkedList* headers = create_linked_list();
	load_to_linkedList(headers, argv[2]);
	printf("Filtering Biostr ... OK\n");
	Biostr_filter(sliced_seqs, headers);
	printf("Write to file ... OK\n");
	Biostr_write_to_file(sliced_seqs, "out.fasta");
	Biostr_free(sliced_seqs);
	printf("-----------------------------------------------------\n");
	struct linkedList* ll = create_linked_list();
	add_node(ll, ">Seq1 skjn species");
	add_node(ll, ">Seq2 dsih species");
	add_node(ll, ">Seq3 dhju species");
	printf("Printing header from linkedList structure\n");
	printf("%s\n", ll->head->header);
	printf("%s\n", ll->head->next->header);
	printf("%s\n", ll->head->next->next->header);
	printf("Using print_linked_list() function\n");
	print_linked_list(ll);
	printf("Inserting at position 1\n");
	insert_at_position(ll, ">Seq4 hdtg species", 1);
	print_linked_list(ll);
	printf("Freeing linkedList*\n");
	free_linked_list(ll);
	printf("----------------------------------------------\n");
	printf("Loading linkedList from file\n");
	struct linkedList* ll2 = create_linked_list();
	load_to_linkedList(ll2, argv[2]);
	print_linked_list(ll2);
	free_linked_list(ll2);
	return 0;
}
