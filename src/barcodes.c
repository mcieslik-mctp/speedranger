#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "barcodes.h"

static int bcinfo_cmp(const void *v1, const void *v2) {
	const BarcodeInfo *b1 = (BarcodeInfo *)v1;
	const BarcodeInfo *b2 = (BarcodeInfo *)v2;
	const bc_t bc1 = b1->bc;
	const bc_t bc2 = b2->bc;
	return (bc1 > bc2) - (bc1 < bc2);
}

void wl_read(BarcodeDict *bcdict, FILE *wl_file)
{
	char buf[1024];

	const size_t n_lines = count_lines(wl_file);
        
	size_t wl_size = 0;

	BarcodeInfo *whitelist = safe_malloc(n_lines * sizeof(*whitelist));

	while (fgets(buf, sizeof(buf), wl_file)) {
		if (strchr(buf, '#'))
			continue;

		whitelist[wl_size].bc = encode_bc(buf);
		whitelist[wl_size].count = 0;
		++wl_size;
	}


	qsort(whitelist, wl_size, sizeof(*whitelist), bcinfo_cmp);

	uint32_t *wl_jumpgate = safe_malloc(POW_2_24 * sizeof(*wl_jumpgate));

	wl_jumpgate[0] = 0;
	uint32_t last_hi = 0;
	for (size_t i = 0; i < wl_size; i++) {
		const bc_t bc = whitelist[i].bc;

		const uint32_t hi = HI24(bc);

		if (hi != last_hi) {
			assert(hi > last_hi);

			for (size_t j = (last_hi + 1); j <= hi; j++)
				wl_jumpgate[j] = i;

			last_hi = hi;
		}
	}

	if (last_hi != 0xFFFFFF) {
		for (size_t j = (last_hi + 1); j < POW_2_24; j++)
			wl_jumpgate[j] = wl_size;
	}

	bcdict->jumpgate = wl_jumpgate;
	bcdict->entries = whitelist;
	bcdict->size = wl_size;
	bcdict->unfound = 0;
}

void wl_dealloc(BarcodeDict *bcdict)
{
	free(bcdict->jumpgate);
	free(bcdict->entries);
}

BarcodeInfo *wl_lookup(BarcodeDict *bcdict, bc_t key)
{
	const uint32_t *wl_jumpgate = bcdict->jumpgate;
	const BarcodeInfo *whitelist = bcdict->entries;
	const size_t wl_size = bcdict->size;

	const uint32_t bc_hi = HI24(key);

	const uint32_t lo = wl_jumpgate[bc_hi];

	if (lo == wl_size) {
		return NULL;
	}

	const uint32_t hi = (bc_hi == 0xFFFFFF ? wl_size : wl_jumpgate[bc_hi + 1]);

	if (lo == hi) {
		return NULL;
	}

	assert(hi > lo);
	BarcodeInfo findme = (BarcodeInfo){ .bc = key };
	BarcodeInfo *target = bsearch(&findme, &whitelist[lo], (hi - lo), sizeof(*whitelist), bcinfo_cmp);
	return target;
}

int wl_increment(BarcodeDict *bcdict, bc_t key)
{
	BarcodeInfo *bcinfo = wl_lookup(bcdict, key);

	if (bcinfo != NULL) {
		++(bcinfo->count);
		return 1;
	} else {
		++(bcdict->unfound);
		return 0;
	}
}

void wl_compute_priors(BarcodeDict *bcdict)
{
	uint64_t total = 0;
	BarcodeInfo *whitelist = bcdict->entries;
	const size_t size = bcdict->size;

	for (size_t i = 0; i < size; i++) {
		total += whitelist[i].count + 1;  // +1 as a pseudocount
	}

	for (size_t i = 0; i < size; i++) {
		whitelist[i].prior = (whitelist[i].count + 1.0)/total;
	}
}

int wl_get_bucket(BarcodeDict *bcdict, BarcodeInfo *bc, const int n_buckets)
{
	return ((bc - bcdict->entries)*n_buckets)/(bcdict->size);
}

void wl_serialize(BarcodeDict *bcdict, FILE *out)
{
	uint32_t *jumpgate = bcdict->jumpgate;
	BarcodeInfo *entries = bcdict->entries;
	const size_t size = bcdict->size;

	for (size_t i = 0; i < POW_2_24; i++) {
		serialize_uint32(out, jumpgate[i]);
	}

	serialize_uint64(out, size);

	for (size_t i = 0; i < size; i++) {
		serialize_uint32(out, entries[i].bc);
		serialize_uint32(out, entries[i].count);
	}
}

void wl_deserialize(BarcodeDict *bcdict, FILE *in)
{
	uint32_t *jumpgate = safe_malloc(POW_2_24 * sizeof(*jumpgate));

	for (size_t i = 0; i < POW_2_24; i++) {
		jumpgate[i] = read_uint32(in);
	}

	const size_t size = read_uint64(in);
	BarcodeInfo *entries = safe_malloc(size * sizeof(*entries));

	for (size_t i = 0; i < size; i++) {
		entries[i].bc = read_uint32(in);
		entries[i].count = read_uint32(in);
	}

	bcdict->jumpgate = jumpgate;
	bcdict->entries = entries;
	bcdict->size = size;
	bcdict->unfound = 0;
}

int correct_barcode(char *barcode, char *barcode_qual, BarcodeDict *wl)
{
#define ILLUMINA_QUAL_OFFSET 33

	uint8_t quals[BC_LEN];

	int n_count = 0;
	for (size_t i = 0; i < BC_LEN; i++) {
		barcode[i] = toupper(barcode[i]);
		quals[i] = barcode_qual[i] - ILLUMINA_QUAL_OFFSET;

		if (!IS_ACGT(barcode[i])) {
			++n_count;
		}
	}

	const bc_t bc0 = ((n_count == 0) ? encode_bc(barcode) : 0);
	BarcodeInfo *bcinfo = ((n_count == 0) ? wl_lookup(wl, bc0) : NULL);

#define CAND_BUF_SIZE (120 * 4 * 4)
	struct bc_string { char bc_str[BC_LEN]; };
	struct bc_string bc_cands[CAND_BUF_SIZE];
	double bc_cand_probs[CAND_BUF_SIZE];
	size_t n_cands = 0;
#undef CAND_BUF_SIZE

	if (bcinfo == NULL) {
		if (n_count > 1) {
			return 0;
		}

		/* examine Hamming-1 neighbors */
		for (size_t i = 0; i < BC_LEN; i++) {
			const char prev = barcode[i];

			if (n_count > 0 && IS_ACGT(prev))
				continue;

			for (size_t j = 0; j < 4; j++) {
				const char new = "ACGT"[j];

				if (new == prev)
					continue;

				barcode[i] = new;
				const bc_t bc = encode_bc(barcode);
				bcinfo = wl_lookup(wl, bc);

				if (bcinfo != NULL) {
					const double prior = bcinfo->prior;
					const double edit_log_prob = MIN(33.0, (double)(quals[i]));
					const double edit_prob = pow(10.0, (-edit_log_prob/10.0));
					const double p = prior*edit_prob;

					memcpy(bc_cands[n_cands].bc_str, barcode, BC_LEN);
					bc_cand_probs[n_cands] = p;
					++n_cands;
				}
			}

			barcode[i] = prev;
		}
	} else {
		assert(n_count == 0);
		memcpy(bc_cands[n_cands].bc_str, barcode, BC_LEN);
		bc_cand_probs[n_cands] = bcinfo->prior;
		++n_cands;

		/* examine Hamming-2 neighbors */
		for (size_t i1 = 0; i1 < BC_LEN; i1++) {
			const char prev1 = barcode[i1];

			for (size_t j1 = 0; j1 < 4; j1++) {
				const char new1 = "ACGT"[j1];

				if (new1 == prev1)
					continue;

				barcode[i1] = new1;

				for (size_t i2 = i1 + 1; i2 < BC_LEN; i2++) {
					const char prev2 = barcode[i2];

					for (size_t j2 = 0; j2 < 4; j2++) {
						const char new2 = "ACGT"[j2];

						if (new2 == prev2)
							continue;

						barcode[i2] = new2;
						const bc_t bc = encode_bc(barcode);
						bcinfo = wl_lookup(wl, bc);

						if (bcinfo != NULL) {
							const double prior = bcinfo->prior;

							const double e1 = MAX(3.0, (double)(quals[i1]) - 1.0);
							const double e2 = MAX(3.0, (double)(quals[i2]) - 1.0);

							const double edit1_log_prob = MIN(33.0, e1);
							const double edit2_log_prob = MIN(33.0, e2);

							const double edit_prob = pow(10.0, (-edit1_log_prob/10.0)) * pow(10.0, (-edit2_log_prob/10.0));
							const double p = prior*edit_prob;

							memcpy(bc_cands[n_cands].bc_str, barcode, BC_LEN);
							bc_cand_probs[n_cands] = p;
							++n_cands;
						}
					}

					barcode[i2] = prev2;
				}

				barcode[i1] = prev1;
			}
		}
	}

	if (n_cands > 0) {
		double total = bc_cand_probs[0];
		size_t max = 0;

		for (size_t i = 1; i < n_cands; i++) {
			total += bc_cand_probs[i];

			if (bc_cand_probs[i] > bc_cand_probs[max])
				max = i;
		}

		if (bc_cand_probs[max]/total > BC_CONF_THRESH) {
			memcpy(barcode, bc_cands[max].bc_str, BC_LEN);
			return 1;
		}
	}

	return 0;

#undef ILLUMINA_QUAL_OFFSET
}
