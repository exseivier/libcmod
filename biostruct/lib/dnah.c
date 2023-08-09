#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "./dnah.h"

void remove_returnc(char* str)
{
	int slen = strlen(str);
	if (str[slen-1] == '\n')
	{
		str[slen-1] = '\0';
	}
}

struct lkdList* loadLines_lkdList(char* filename)
{
	//	Assert that header sizes less than 255 characters.
	//int buffer = BUFFER;
	char carrier[256];
	FILE* f = fopen(filename, "r");
	struct lkdList* head = (struct lkdList*)calloc(sizeof(struct lkdList), 1);
	struct lkdList* tail;
	fgets(carrier, 256, f);
	remove_returnc(carrier);
	head->header = strdup(carrier);
	head->next = NULL;
	tail = head;
	while (fgets(carrier, 256, f))
	{
		struct lkdList* node = (struct lkdList*)malloc(sizeof(struct lkdList));
		remove_returnc(carrier);
		node->header = strdup(carrier);
		node->next = NULL;
		tail->next = node;
		tail = node;
	}
	fclose(f);
	return head;
}

struct DNA* loadDNASeqs(char* filename)
{
	long int fileSize = 0;
	FILE* fh = fopen(filename, "rb");
	fseek(fh, 0, SEEK_END);
	fileSize = ftell(fh);
	rewind(fh);
	char* fileStr = (char*)calloc(sizeof(char), fileSize);
	fread(fileStr, sizeof(char), fileSize, fh);
	fclose(fh);
	fileStr[fileSize] = '\0';


	int ITEM_BUFFER = 50;
	int SEQ_BUFFER = 1048576;

	int stCount = 0;
	//int endCount = 0;
	int itemCount = -1;
	int lineSize = 0;
	int seq_len = 0;
	struct DNA* items = (struct DNA*)calloc(sizeof(struct DNA), 1);
	items->ids = (char**)calloc(sizeof(char*), ITEM_BUFFER);
	items->start = (int*)calloc(sizeof(int), ITEM_BUFFER);
	items->end = (int*)calloc(sizeof(int), ITEM_BUFFER);
	items->seq_len = (int*)calloc(sizeof(int), ITEM_BUFFER);
	items->seq = (char*)calloc(sizeof(char), SEQ_BUFFER);
		
	
	char* line = strtok(fileStr, "\n");
	while (line != NULL)
	{
		if (line[0] == '>')
		{
			if (seq_len > 0)
			{
				items->end[itemCount] = seq_len-1;
			}
			itemCount++;
			if (itemCount == ITEM_BUFFER)
			{
				ITEM_BUFFER = ITEM_BUFFER * 2;
				items->ids = (char**)realloc(items->ids, sizeof(char*) * ITEM_BUFFER);
				items->start = (int*)realloc(items->start, sizeof(int) * ITEM_BUFFER);
				items->end = (int*)realloc(items->end, sizeof(int) * ITEM_BUFFER);
				items->seq_len = (int*)realloc(items->seq_len, sizeof(int) * ITEM_BUFFER);
			}
			lineSize = strlen(line);
			items->ids[itemCount] = (char*)calloc(sizeof(char), lineSize);
			items->ids[itemCount] = strdup(line);
			stCount = seq_len;
			items->start[itemCount] = stCount;
		}
		else
		{
			lineSize = strlen(line);
			if(seq_len + lineSize >= SEQ_BUFFER)
			{
				SEQ_BUFFER = SEQ_BUFFER * 2;
				items->seq = (char*)realloc(items->seq, sizeof(char) * SEQ_BUFFER);
			}


			int j = 0;
			for (int i = seq_len; i < seq_len + lineSize; i++)
			{
				items->seq[i] = line[j];
				j++;
			}
			seq_len = seq_len + lineSize;
		}
		line = strtok(NULL, "\n");
	}
	items->end[itemCount] = seq_len-1;
	itemCount++;
	items->ids[itemCount] = NULL;
	int i = 0;
	items->hide = (int*)calloc(sizeof(int), ITEM_BUFFER);
	while (items->ids[i] != NULL)
	{
		items->hide[i] = 0;
		items->seq_len[i] = items->end[i] - items->start[i] + 1;
		i++;
	}
	return items;
}

void freeDNA(struct DNA* items)
{
	free(items->seq);
	free(items->start);
	int i = 0;
	while (items->ids[i] != NULL)
	{
		free(items->ids[i]);
		i++;
	}
}

void freeLkdList(struct lkdList* lkdlst)
{
	struct lkdList* tmp = lkdlst;
	struct lkdList* aux;
	while (tmp != NULL)
	{
		free(tmp->header);
		aux = tmp;
		tmp = tmp->next;
		free(aux);
	}
}

struct DNA* splitBioString(struct DNA* seqs, int size, int step)
{
	int ITEMS_BUFFER = 0;
	int itemCount = 0;
	int seq_len = strlen(seqs->seq);
	struct DNA* sptSeqs = (struct DNA*)calloc(sizeof(struct DNA), 1);
	sptSeqs->seq = (char*)calloc(sizeof(char), seq_len+1);
	for (int i = 0; i < seq_len; i++)
	{
		sptSeqs->seq[i] = seqs->seq[i];
	}
	sptSeqs->seq[seq_len] = '\0';

	int i = 0;
	char header_buffer[256];
	double next_items = 0.0;
	while (seqs->ids[i] != NULL)
	{
		next_items = (double)(seqs->end[i] - seqs->start[i] + 1) / step;
		ITEMS_BUFFER = (int)ceil(itemCount + next_items);
		sptSeqs->ids = (char**)realloc(sptSeqs->ids, sizeof(char*) * ITEMS_BUFFER);
		sptSeqs->start = (int*)realloc(sptSeqs->start, sizeof(int) * ITEMS_BUFFER);
		sptSeqs->end = (int*)realloc(sptSeqs->end, sizeof(int) * ITEMS_BUFFER);
		sptSeqs->seq_len = (int*)realloc(sptSeqs->seq_len, sizeof(int) * ITEMS_BUFFER);
		sptSeqs->hide = (int*)realloc(sptSeqs->hide, sizeof(int) * ITEMS_BUFFER);
		
		for(int j = seqs->start[i]; j < seqs->end[i]; j = j+step)
		{
			sprintf(header_buffer, "%s_%d", seqs->ids[i], itemCount);
			sptSeqs->ids[itemCount] = strdup(header_buffer);
			sptSeqs->start[itemCount] = j;
			sptSeqs->end[itemCount] = (seqs->end[i] < (j+size-1)) ? seqs->end[i] : (j+size-1);
			sptSeqs->seq_len[itemCount] = sptSeqs->end[itemCount] - sptSeqs->start[itemCount] + 1;
			sptSeqs->hide[itemCount] = 0;
			itemCount++;
		}
		i++;
	}

	if (itemCount > ITEMS_BUFFER)
	{
		printf("Can't be possible. That means a memory allocation error!");
		exit(1);
	}
	if (itemCount == ITEMS_BUFFER)
	{
		ITEMS_BUFFER = ITEMS_BUFFER + 1;
		sptSeqs->ids = (char**)realloc(sptSeqs->ids, sizeof(char*) * ITEMS_BUFFER);
		sptSeqs->start = (int*)realloc(sptSeqs->start, sizeof(int) * ITEMS_BUFFER);
		sptSeqs->end = (int*)realloc(sptSeqs->end, sizeof(int) * ITEMS_BUFFER);
		sptSeqs->seq_len = (int*)realloc(sptSeqs->seq_len, sizeof(int) * ITEMS_BUFFER);
		sptSeqs->hide = (int*)realloc(sptSeqs->hide, sizeof(int) * ITEMS_BUFFER);
	}
	sptSeqs->ids[itemCount] = '\0';
	sptSeqs->start[itemCount] = '\0';
	sptSeqs->end[itemCount] = '\0';
	sptSeqs->seq_len[itemCount] = '\0';
	sptSeqs->hide[itemCount] = '\0';

	return sptSeqs;
}

void filterBioseq(struct DNA* bs, struct lkdList* headers)
{
	int itemCount = 0;
	struct lkdList* previous = NULL;
	struct lkdList* current = NULL;
	while (bs->ids[itemCount] != NULL)
	{
		current = headers;
		previous = NULL;
		if (bs->hide[itemCount] == 1)
		{
			itemCount++;
			continue;
		}
		while(current != NULL)
		{
			if (strcmp(bs->ids[itemCount], current->header) == 0)
			{
				printf("%s - %s\n", bs->ids[itemCount], current->header);
				printf("%s\n%d\n", bs->ids[itemCount], bs->seq_len[itemCount]);
				bs->hide[itemCount] = 1;
				if (previous == NULL)
				{
					current = current->next;
					headers = current;
				}
				else
				{
					previous->next = current->next;
					free(current->header);
					current = NULL;
					current = previous->next;
				}
				break;
			}
			else
			{
				//printf("Current header: %s\n", current->header);
				previous = current;
				current = current->next;
			}
		}
		//printf("Next Bioseq: %s\n", bs[bs_count]->id);
		itemCount++;
	}
}

void writeNoHideToFile(struct DNA* bs)
{
	FILE* output = fopen("output.fas", "w");
	int i = 0;
	while (bs->ids[i] != NULL)
	{
		if (bs->hide[i] == 0)
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
			fprintf(output, "%s\n%s\n", bs->ids[i], seqCarrier);
			free(seqCarrier);
			seqCarrier = NULL;
		}
		i++;
	}
	fclose(output);
}

void sampleSeqs(struct DNA* bs, int perc, char* outfile)
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


