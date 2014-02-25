#include "common.h"
/* Header file for utils.c */
#ifdef HAVE_CONFIG_H
#include <config.h>
#ifndef HAVE_STRDUP
CPPEXTERN char *strdup(const char *s);
#endif
#endif
#ifdef WITH_DMALLOC
/* use dmalloc library to check for memory management bugs */
#include "dmalloc.h"
#define space(S) calloc(1,(S))
#else
CPPEXTERN /*@only@*/ /*@notnull@*/
void  *space(unsigned size) /*@ensures MaxSet(result) == (size-1);@*/;
			    /* allocate space safely */
CPPEXTERN /*@only@*/ /*@notnull@*/
void  *xrealloc(/*@null@*/ /*@only@*/ /*@out@*/ /*@returned@*/ void *p, unsigned size) /*@modifies *p @*/ /*@ensures MaxSet(result) == (size-1) @*/;
#endif

CPPEXTERN /*@exits@*/ void nrerror(const char message[]);  /* die with error message */
CPPEXTERN void   init_rand(void);                /* make random number seeds */
CPPEXTERN unsigned short xsubi[3];               /* current 48bit random number */
CPPEXTERN double urn(void);                      /* random number from [0..1] */
CPPEXTERN int    int_urn(int from, int to);      /* random integer */
CPPEXTERN void   filecopy(FILE *from, FILE *to); /* inefficient `cp' */
CPPEXTERN /*@observer@*/ char  *time_stamp(void);               /* current date in a string */
CPPEXTERN /*@only@*/ /*@notnull@*/ char  *random_string(int l, const char symbols[]);
/* random string of length l using characters from symbols[] */
CPPEXTERN int    hamming(const char *s1, const char *s2);
/* calculate hamming distance */
CPPEXTERN /*@only@*/ /*@null@*/ char  *get_line(const FILE *fp); /* read one (arbitrary length) line from fp */


CPPEXTERN char *pack_structure(const char *struc);
/* pack secondary secondary structure, 5:1 compression using base 3 encoding */
CPPEXTERN char *unpack_structure(const char *packed);
/* unpack sec structure packed with pack_structure() */
CPPEXTERN short *make_pair_table(const char *structure);
/* returns a newly allocated table, such that:  table[i]=j if (i.j) pair or
   0 if i is unpaired, table[0] contains the length of the structure. */

CPPEXTERN int bp_distance(const char *str1, const char *str2);
/* dist = {number of base pairs in one structure but not in the other}
   same as edit distance with open-pair close-pair as move-set */
