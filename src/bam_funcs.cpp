#include <Rcpp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sam.h>

using namespace Rcpp;

//                 barcode,     chromosome,  beg,          end,          n
typedef std::tuple<std::string, std::string, std::int32_t, std::int32_t, std::int32_t> Molecule;


//' Extract putative DNA molecules from barcoded BAM file.
//'
//' @param bam_fn filename
//' @param tagname barcode tag name
//' @param max_dist maximum distance between reads in a molecule
//' @export
// [[Rcpp::export]]
DataFrame extractMolecules(std::string bam_fn, std::string tagname="BX", std::int32_t max_dist=50e4) {
  const char *tagname_c = tagname.c_str();
  std::vector<Molecule> complete_mols;
  std::unordered_map<std::string,Molecule> current_mols;
  std::uint8_t *tagid;
  std::string barcode;
  std::string chromosome;
  std::int32_t beg;
  std::int32_t end;
  Molecule molecule;

  
  Rcpp::Rcout << "Started reading bam file" << std::endl;
  // open bam file
  samFile *fp_in = hts_open(bam_fn.c_str(), "r");
  // header parse
  bam_hdr_t *bamHdr = sam_hdr_read(fp_in);
  // initialize an alignment
  bam1_t *rec = bam_init1(); 
  while(sam_read1(fp_in, bamHdr, rec) > 0) {
    // skip some reads
    if ( rec->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
    tagid = bam_aux_get(rec, tagname_c);
    if (tagid) {
      barcode = bam_aux2Z(tagid);
      chromosome = bamHdr->target_name[rec->core.tid];
      beg = rec->core.pos + 1;
      end = bam_endpos(rec);
      if (current_mols.count(barcode) != 0) {
        // molecule found
        molecule = current_mols[barcode];
        if ( (std::get<1>(molecule) == chromosome) && (beg - std::get<3>(molecule) < max_dist) ) {
          // extend molecule
          std::get<3>(molecule) = std::max(end, std::get<3>(molecule));
          std::get<4>(molecule)++;
        } else {
          // complete molecule
          complete_mols.push_back(molecule);
          // create new molecule
          molecule = std::make_tuple(barcode, chromosome, beg, end, 1);
        }
      } else {
        // molecule not found
        // create new molecule
        molecule = std::make_tuple(barcode, chromosome, beg, end, 1);
      }
      current_mols[barcode] = molecule;
    }
  }
  bam_destroy1(rec);
  bam_hdr_destroy(bamHdr);
  sam_close(fp_in);
  Rcpp::Rcout << "Finished reading bam file" << std::endl;
  // complete all remaining molecules
  Rcpp::Rcout << "Started molecule completion" << std::endl;
  for(std::unordered_map<std::string,Molecule>::iterator it = current_mols.begin(); it != current_mols.end();) {
    complete_mols.push_back(it->second);
    it = current_mols.erase(it);
  }
  Rcpp::Rcout << "Finished molecule completion" << std::endl;
  // pass results to R
  Rcpp::Rcout << "Started data copy" << std::endl;
  std::vector<std::string> res_barcode;
  std::vector<std::string> res_chromosome;
  std::vector<std::int32_t> res_beg;
  std::vector<std::int32_t> res_end;
  std::vector<std::int32_t> res_n;
  std::vector<Molecule>::iterator it = complete_mols.end();
  while(!complete_mols.empty()) {
    molecule = complete_mols.back();
    res_barcode.push_back(std::get<0>(molecule));
    res_chromosome.push_back(std::get<1>(molecule));
    res_beg.push_back(std::get<2>(molecule));
    res_end.push_back(std::get<3>(molecule));
    res_n.push_back(std::get<4>(molecule));
    complete_mols.pop_back();
  }
  Rcpp::Rcout << "Finished data copy" << std::endl;
  DataFrame res = DataFrame::create(
                               Named("chr")=res_chromosome,
                               Named("beg")=res_beg,
                               Named("end")=res_end,
                               Named("barcode")=res_barcode,
                               Named("n")=res_n
                               );
  return res;
}

// template<std::size_t> struct int_{};

// template <class Tuple, size_t Pos>
// std::ostream& print_tuple(std::ostream& out, const Tuple& t, int_<Pos> ) {
//   out << std::get< std::tuple_size<Tuple>::value-Pos >(t) << ',';
//   return print_tuple(out, t, int_<Pos-1>());
// }

// template <class Tuple>
// std::ostream& print_tuple(std::ostream& out, const Tuple& t, int_<1> ) {
//   return out << std::get<std::tuple_size<Tuple>::value-1>(t);
// }

// template <class... Args>
// std::ostream& operator<<(std::ostream& out, const std::tuple<Args...>& t) {
//   out << '(';
//   print_tuple(out, t, int_<sizeof...(Args)>());
//   return out << ')';
// }
