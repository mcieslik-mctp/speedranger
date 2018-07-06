#ifndef BARCODES_H
#define BARCODES_H

#include "util.h"

/*
 * Useful structures for computing barcode information
 */

typedef struct {
	bc_t bc;
	uint32_t count;
	double prior;
} BarcodeInfo;

typedef struct {
	uint32_t *jumpgate;
	BarcodeInfo *entries;
	size_t size;
	uint32_t unfound;
} BarcodeDict;

/* misc. bc functions */
void wl_read(BarcodeDict *bcdict, FILE *wl_file);
void wl_dealloc(BarcodeDict *bcdict);
BarcodeInfo *wl_lookup(BarcodeDict *bcdict, bc_t key);
int wl_increment(BarcodeDict *bcdict, bc_t key);
void wl_compute_priors(BarcodeDict *bcdict);
int wl_get_bucket(BarcodeDict *bcdict, BarcodeInfo *bc, const int n_buckets);

/* serialize/deserializes barcode dictionaries */
void wl_serialize(BarcodeDict *bcdict, FILE *out);
void wl_deserialize(BarcodeDict *bcdict, FILE *in);

/* barcode correction */
int correct_barcode(char *barcode, char *barcode_qual, BarcodeDict *wl);

/* corrects barcodes and generate new FASTQs */
/* This is essentially a translation of Long Ranger's barcode correction scheme. */
void preprocess_fastqs(const char *cts, const char *inp, const char *out);

#define BC_CONF_THRESH	0.975
#define BC_LEN		16
#define MATE1_TRIM	7
#define MIN_READ_LEN	15

#endif /* BARCODES_H */
