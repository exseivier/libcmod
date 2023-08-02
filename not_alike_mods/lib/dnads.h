/*
 * DNA data structures C header
 * Author: Javier Montalvo-Arredondo.
 * Contact: buitrejma@gmail.com
 * UNIVERSIDAD AUTONOMA AGRARIA ANTONIO NARRO
 * DEPARTAMENTO DE CIENCIAS BÁSICAS
 * [BASIC SCIENCES DEPARTMENT]
 */

#ifndef _DNADS_H_
#define _DNADS_H_

#define BUFFER 2048
#define SEQ_LEN_BUFFER 1048576

//	Bioseq data structure
//	Stores the id, sequence and hide information.
struct Bioseq
{
	char* id;
	char* seq;
	int hide;
};

//	Linked list data structure used to hold ids.
struct lkdList
{
	char* header;
	struct lkdList* next;	
};

//	Computes the number of sequences a Bioseq structure has.
int howManySeqs(struct Bioseq** bs);

//	Removes the last return character ('\n') of a string.
void remove_returnc(char* str);

//	Loads the sequences from FASTA file and stores them into a Bioseq structure.
struct Bioseq** loadBioSeq(char* filename);

//	Free dynamically allocated memory used to create a Bioseq structure.
void freeBioseq(struct Bioseq** bs);

//	Extracts a subsequence from a string.
char* extractSubSeq(char* str, int start, int end);

//	Split sequences stored in Bioseq object and the resulting
//	split sequences are stored in another Bioseq object.
struct Bioseq** splitBioseq(struct Bioseq** bs, int size, int step);


//	Free dynamically allocated memory used to create a linked list of sequences ids.
void freeLkdList(struct lkdList* lkdlst);

//	Creates a linked list of sequences id from a header.txt file.
/*
 * 	Example of header.txt file.
 * 	---------------------------
 * 	>Seq_1
 * 	>Seq_2
 * 	>Seq_3
 * 	...
 * 	>Seq_N
 * 	EOF
 */
struct lkdList* loadLines_lkdList(char* filename);


//	Filters elements in Bioseq structure.
//	If the id is present in header.txt file, hide information is set to true (1).
void filterBioseq(struct Bioseq** bs, struct lkdList* headers);


//	Writes to file the sequences stored in Bioseq structure only if hide information is set to False (0).
void writeNoHideToFile(struct Bioseq** bs);

#endif
