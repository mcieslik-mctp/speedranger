#include <Rcpp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sam.h>
#include <zlib.h>
#include <stdint.h>
using namespace Rcpp;


extern "C" 
{
  
  #include "barcodes.h"
  #include "seqtk/kseq.h"
  KSEQ_INIT(gzFile, gzread)
  
}


//' Count barcodes
//'
//' @param inp_fn input FQ file
//' @param out_fn output counts file
//' @param wl_fn output counts file
//' @export
// [[Rcpp::export]]
StringVector countBarcodes(std::string inp_fn, std::string out_fn, std::string wl_fn) {

  if (access(inp_fn.c_str(), F_OK) == -1) {
    stop("inp_fn could not be opened for reading");
  }

  FILE *out_file = fopen(out_fn.c_str(), "wb");
  if (out_file == NULL) {
    stop("out_fn could not be opened for writing");
  }

  if (access(wl_fn.c_str(), F_OK) == -1) {
    stop("wl_fn could not be opened for reading");
  }
  
  FILE *wl_file = fopen(wl_fn.c_str(), "rb");
  BarcodeDict wldict;
  wl_read(&wldict, wl_file);

  // open file handles
  gzFile inp_file;
  inp_file = gzopen(inp_fn.c_str(), "rb");
  kseq_t *ks_inp;
  ks_inp = kseq_init(inp_file);

  char barcode[BC_LEN+1] = {'\0'};
  while (kseq_read(ks_inp) >= 0) {
    for (size_t i = 0; i < BC_LEN; i++) {
      if (IS_ACGT(ks_inp->seq.s[i])) {
        barcode[i] = ks_inp->seq.s[i];
      } else {
        barcode[0] = '\0';
        break;
      }
    }
    if (barcode[0] != '\0') {
      const bc_t bc = encode_bc(barcode);
      wl_increment(&wldict, bc);
    }
  }
  // count_barcodes(&wldict, inp_file);
  wl_serialize(&wldict, out_file);
  /* clean-up */
  kseq_destroy(ks_inp);
  wl_dealloc(&wldict);
  gzclose(inp_file);
  fclose(out_file);
  fclose(wl_file);
  return out_fn;
}


//' Preprocess barcoded FASTQ files.
//'
//' @param cts_inp_fn input FQ file read 1
//' @param fq1_inp_fn input FQ file read 1
//' @param fq2_inp_fn input FQ file read 2
//' @param fq1_out_fn output FQ file read 1
//' @param fq2_out_fn output FQ file read 2
//' @export
// [[Rcpp::export]]
StringVector preprocessFastq(std::string cts_inp_fn, std::string fq1_inp_fn, std::string fq2_inp_fn,
                                                     std::string fq1_out_fn, std::string fq2_out_fn) {
  
  if (access(cts_inp_fn.c_str(), F_OK) == -1) {
    stop("cts_inp_fn could not be opened for reading");
  }
  
  if(
     (access(fq1_inp_fn.c_str(), F_OK) == -1) ||
     (access(fq2_inp_fn.c_str(), F_OK) == -1)
     ) {
    stop("fq1_inp_fn or fq2_inp_fn could not be opened for reading");
  }

  FILE *fq_out[2];
  fq_out[0] = fopen(fq1_out_fn.c_str(), "wb");
  fq_out[1] = fopen(fq2_out_fn.c_str(), "wb");

  if (fq_out[0] == NULL || fq_out[1] == NULL) {
    stop("fq1_out_fn or fq2_out_fn could not be opened for writing");
  }

  // read counts
  FILE *cts_file = fopen(cts_inp_fn.c_str(), "rb");
  BarcodeDict wl;
  wl_deserialize(&wl, cts_file);
  wl_compute_priors(&wl);
  
  // open file handles
  gzFile fq_inp[2];
  fq_inp[0] = gzopen(fq1_inp_fn.c_str(), "rb");
  fq_inp[1] = gzopen(fq2_inp_fn.c_str(), "rb");

  kseq_t *ks_inp[2];
  ks_inp[0] = kseq_init(fq_inp[0]);
  ks_inp[1] = kseq_init(fq_inp[1]);

  kseq_t *r1, *r2;
  FILE *f1, *f2;

  char barcode[BC_LEN+1] = {'\0'};
  char barcode_qual[BC_LEN+1] = {'\0'};
  while ((kseq_read(ks_inp[0]) >= 0) && (kseq_read(ks_inp[1]) >= 0)) {

    f1 = fq_out[0];
    r1 = ks_inp[0];
    f2 = fq_out[1];
    r2 = ks_inp[1];

    if (r1->seq.l < (BC_LEN + MATE1_TRIM + MIN_READ_LEN)) {
      continue;
    }

    if (r2->seq.l < (MIN_READ_LEN)) {
      continue;
    }
    
    memcpy(barcode, r1->seq.s, BC_LEN);
    memcpy(barcode_qual, r1->qual.s, BC_LEN);
    const int good_barcode = correct_barcode(barcode, barcode_qual, &wl);

    // read1
    fputc('@', f1);
    fputs(r1->name.s, f1);
    fputs("\tBX:Z:", f1); fputs(barcode, f1); fputs("-1", f1);
    fputs("\tCX:Z:", f1); fputc(good_barcode ? '1' : '0', f1);
    fputc('\n', f1); fputs(r1->seq.s+BC_LEN+MATE1_TRIM, f1);
    fputs("\n+\n", f1); fputs(r1->qual.s+BC_LEN+MATE1_TRIM, f1);
    fputc('\n', f1);

    // read2
    fputc('@', f2);
    fputs(r2->name.s, f2);
    fputs("\tBX:Z:", f2); fputs(barcode, f2); fputs("-1", f2);
    fputs("\tCX:Z:", f2); fputc(good_barcode ? '1' : '0', f2);
    fputc('\n', f2); fputs(r2->seq.s, f2);
    fputs("\n+\n", f2); fputs(r2->qual.s, f2);
    fputc('\n', f2);
      
  }
  
  // cleanup
  fclose(cts_file);  
  kseq_destroy(ks_inp[0]);
  kseq_destroy(ks_inp[1]);
  gzclose(fq_inp[0]);
  gzclose(fq_inp[1]);
  fclose(fq_out[0]);
  fclose(fq_out[1]);

  // output
  StringVector res(2);
  res[0] = fq1_out_fn;
  res[1] = fq2_out_fn;
  return res;
}
