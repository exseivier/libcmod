#ifndef _DNAH_H_
#define _DNAH_H_

/*
 * DNA data structures C header.
 * Author: Javier Montalvo-Arredondo.
 * Contact: buitrejma@gmail.com
 * Universidad Autonoma Agraria Antonio Narro.
 * Departamento de Ciencias Basicas.
 * Basic Sciences Department.
 */


struct DNA
{
	char* seq;
	char** ids;
	int* start;
	int* end;
	int* seq_len;
	int* hide;
};

struct lkdList
{
	char* header;
	struct lkdList* next;	
};

//	Loads DNA sequences from FASTA file and stores into DNA structure.
struct DNA* loadDNASeqs(char* filename);

//	Free dynamically allocated memory used to hold a DNA structure.
void freeDNA(struct DNA* items);

//	Free dynamically allocated memory used to hold a linked list structure.
void freeLkdList(struct lkdList* lkdlst);

//	Removes last \n character from reading file buffer.
void remove_returnc(char* str);

//	Creates a linked list with headers of sequences. 
struct lkdList* loadLines_lkdList(char* filename);

//	Splits DNA string in subsequences of determined size at defined step.
struct DNA* splitBioString(struct DNA* seqs, int size, int step);

//	Filters subsequences which has a BLAST hit.
void filterBioseq(struct DNA* bs, struct lkdList* headers);

//	Writes to file those sequences which has not BLAST hit.
void writeNoHideToFile(struct DNA* bs);

#endif
