#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "./dnads.h"

#define BUFFER 2048
#define SEQ_LEN_BUFFER 1048576

int howManySeqs(struct Bioseq** bs)
{
	int i = 0;
	while (bs[i] != NULL)
	{
		i++;
	}
	return i;
}

void remove_returnc(char* str)
{
	int slen = strlen(str);
	if (str[slen-1] == '\n')
	{
		str[slen-1] = '\0';
	}
}

struct Bioseq** loadBioSeq(char* filename)
{
	//	Load file content to an array.
	
	int BioseqSize = 50;
	FILE* input = fopen(filename, "rb");
	if (NULL == input)
	{
		printf("I couldn't open the file...\n");
		exit(1);
	}

	fseek(input, 0, SEEK_END);
	long int fileSize = ftell(input);
	rewind(input);

	char* strArray = (char*)calloc(sizeof(char), fileSize);
	fread(strArray, sizeof(char), fileSize, input);	
	fclose(input);
	strArray[fileSize] = '\0';

	//	Store the sequences into Bioseq.
	struct Bioseq** seqs = (struct Bioseq**)calloc(sizeof(struct Bioseq*), BioseqSize);
	int num_s = -1;
	int seq_len = 0;
	int seq_len_buffer = SEQ_LEN_BUFFER;
	char* line = strtok(strArray, "\n");
	while (NULL != line)
	{
		if (num_s == BioseqSize-1)
		{
			BioseqSize = BioseqSize * 2;
			seqs = realloc(seqs, sizeof(struct Bioseq*) * BioseqSize);
		}
		if (line[0] == '>')
		{
			num_s++;
			seqs[num_s] = calloc(sizeof(struct Bioseq), 1);
			seqs[num_s]->id = strdup(line);
			seqs[num_s]->seq = calloc(sizeof(char), seq_len_buffer);
		}
		else
		{
			seq_len = seq_len + strlen(line);
			if (seq_len >= seq_len_buffer-1)
			{
				seq_len_buffer = seq_len_buffer * 2;
				seqs[num_s]->seq = realloc(seqs[num_s]->seq, sizeof(char) * seq_len_buffer);
			}
			seqs[num_s]->seq = strcat(seqs[num_s]->seq, line);
		}
		line = strtok(NULL, "\n");
	}
	if (num_s == BioseqSize-1)
	{
		BioseqSize = BioseqSize + 1;
		seqs = realloc(seqs, sizeof(struct Bioseq*)*BioseqSize);
	}
	seqs[++num_s] = '\0';
	return seqs;

}

void freeBioseq(struct Bioseq** bs)
{
	int i = 0;
	while (bs[i] != NULL)
	{
		free(bs[i]->id);
		free(bs[i]->seq);
		free(bs[i]);
		i++;
	}
	free(bs);
}

char* extractSubSeq(char* str, int start, int end)
{
	if (end <= start)
	{
		printf("[ERROR] : End position is bigger than start position.\n");
		printf("[TIP] : Make sure that start position be bigger than end position.\n");
		exit(1);
	}
	int size = end-start;
	char* subStr = (char*)calloc(sizeof(char), (size+1));
	//subStr = NULL;  // Don't set NULL to a char*. It causes segmentation fault.
	for(int i = start; i < end; i++)
	{
		subStr[i-start] = str[i];
	}
	subStr[size] = '\0';
	return subStr;
}

struct Bioseq** splitBioseq(struct Bioseq** bs, int size, int step)
{
	int new_buffer = BUFFER;
	char header_buffer[BUFFER];
	struct Bioseq** sptSeqs = (struct Bioseq**)calloc(sizeof(struct Bioseq*), new_buffer);
	// struct Bioseq** bs counter.
	int i = 0;
	//	Sequence length.
	int slen = 0;
	//	Sequence start.
	int start = 0;
	//	Sequence end.
	int end = size;
	//	Subsequence counter
	int subi = 0;

	while (bs[i] != NULL)
	{
		slen = strlen(bs[i]->seq);
		start = 0;
		end = size;
		while (start <= slen)
		{
			if (subi == new_buffer)
			{
				new_buffer = new_buffer * 2;
				struct Bioseq** tmp = (struct Bioseq**)realloc(sptSeqs, sizeof(struct Bioseq*) * new_buffer);
				if (tmp != NULL) sptSeqs = tmp;
				else
				{
					printf("I can't allocate that amount of memory!...\n");
					freeBioseq(sptSeqs);
					exit(1);
				}
			}
			sprintf(header_buffer, "%s_%d", bs[i]->id, subi);
			sptSeqs[subi] = (struct Bioseq*)calloc(sizeof(struct Bioseq), 1);
			sptSeqs[subi]->id = strdup(header_buffer);
			sptSeqs[subi]->hide = 0;
			sptSeqs[subi]->seq = extractSubSeq(bs[i]->seq, start, end);
			start = start + step;
			end = end + step;
			subi++;
		}
		i++;
	}
	sptSeqs[subi] = NULL;
	return sptSeqs;
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

struct lkdList* loadLines_lkdList(char* filename)
{
	//	Assert that header sizes less than 255 characters.
	//int buffer = BUFFER;
	char carrier[BUFFER];
	FILE* f = fopen(filename, "r");
	struct lkdList* head = (struct lkdList*)malloc(sizeof(struct lkdList));
	struct lkdList* tail;
	fgets(carrier, BUFFER, f);
	remove_returnc(carrier);
	head->header = strdup(carrier);
	head->next = NULL;
	tail = head;
	while (fgets(carrier, BUFFER, f))
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

void filterBioseq(struct Bioseq** bs, struct lkdList* headers)
{
	int bs_count = 0;
	struct lkdList* previous = NULL;
	struct lkdList* current = NULL;
	while (bs[bs_count] != NULL)
	{
		current = headers;
		previous = NULL;
		if (bs[bs_count]->hide == 1)
		{
			bs_count++;
			continue;
		}
		while(current != NULL)
		{
			if (strcmp(bs[bs_count]->id, current->header) == 0)
			{
				printf("%s - %s\n", bs[bs_count]->id, current->header);
				printf("%s\n%s\n", bs[bs_count]->id, bs[bs_count]->seq);
				bs[bs_count]->hide = 1;
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
		bs_count++;
	}
}

void writeNoHideToFile(struct Bioseq** bs)
{
	FILE* output = fopen("output.fas", "w");
	int i = 0;
	while (bs[i] != NULL)
	{
		if (bs[i]->hide == 0)
		{
			fprintf(output, "%s\n%s\n", bs[i]->id, bs[i]->seq);
		}
		i++;
	}
	fclose(output);
}


