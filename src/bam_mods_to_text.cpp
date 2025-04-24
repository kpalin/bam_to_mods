/*

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

// #include <config.h>

#include <ctype.h>
#include <math.h>
#include <stdio.h>

#include <getopt.h>
#include <htslib/faidx.h>
#include <htslib/sam.h>

#include <map>

#include <string>

#include <iostream>

#include <sstream>

#include "version.h"
#include <algorithm>
#include <cassert>
#include <set>
#include <string>
#include <vector>

static int header_flag = 1;
static int split_strand_flag = 0;
static bool combine_hemi = true; // Combine hemi-methylation

static int phase_flag = 1;
static int min_mod_prob = 209; // Ln odds 1.5  >81%
static int max_mod_prob = 46;  // < 18.2%
static int min_mapq = 20;
static int min_baseq = 7;
static int extra_hts_threads = 2;
static int exclude_filter =
    (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY);
// samtools view -F 1796 -u
// /mnt/cgnano/projects/promethion/kpalin/dev/Fam_c461_1_19_0711NK_meg/phase/longshot.Fam_c461_1_19_0711NK_meg.phased.cram
// chr9:170012-170012  |samtools mpileup -Q 7 -q 20 --output-extra PS - |grep
// 170012
typedef struct {
  samFile *fp;
  sam_hdr_t *fp_hdr;
  hts_itr_t *itr;
} plp_dat;

enum WhichStrandOutput {
  STRAND_FORWARD = 0,
  STRAND_REVERSE = 1,
  STRAND_BOTH = -1
};

static int readaln(void *data, bam1_t *b) {
  plp_dat *g = (plp_dat *)data;
  int ret;

  while (1) {
    if (g->itr) {
      ret = sam_itr_next(g->fp, g->itr, b);
    } else {
      ret = sam_read1(g->fp, g->fp_hdr, b);
    }

    if (ret < 0)
      break;
    if (b->core.qual < min_mapq)
      continue;
    if (b->core.flag & exclude_filter)
      continue;
    break;
  }

  return ret;
}

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif


// Initialise and destroy the base modifier state data. This is called
// as each new read is added or removed from the pileups.
int pileup_cd_create(void *data, const bam1_t *b, bam_pileup_cd *cd) {
  hts_base_mod_state *m = hts_base_mod_state_alloc();
  bam_parse_basemod(b, m);
  cd->p = m;
  return 0;
}

int pileup_cd_destroy(void *data, const bam1_t *b, bam_pileup_cd *cd) {
  hts_base_mod_state_free((hts_base_mod_state *)cd->p);
  return 0;
}

#include <regex>

const char complement(char const n) {
  switch (n) {
  case 'A':
    return 'T';
  case 'T':
    return 'A';
  case 'G':
    return 'C';
  case 'C':
    return 'G';
  case 'N':
    return 'N';
  }
  std::cerr << n << '\n';
  assert(false);
  return ' ';
}

std::string complement(std::string &w) {

  std::stringstream c;
  for (std::string::iterator itr = w.begin(); itr != w.end(); ++itr) {
    c << complement((char)(*itr));
  }
  return c.str();
}
std::string revComplement(std::string &fwd) {
  std::stringstream rev;
  for (std::reverse_iterator<std::string::iterator> itr = fwd.rbegin();
       itr != fwd.rend(); ++itr) {
    rev << complement((char)(*itr));
  }
  return rev.str();
}
struct mod_id_code {
private:
  char canonical; // Read original strand base for the modification, e.g. 'C' or
  // 'T'
  int value; // Modification code, e.g. 'm' or 'a' or 123 (As defined in SAM
             // spec)
  bool rev_strand; // Strand of the modification with respect to the original
                   // read strand, e.g. 0 for '+' strand and 1 for '-' strand

public:
  explicit mod_id_code(char canonical_base, int modified_base = 0,
                       bool is_rev = false)
      : canonical(canonical_base), value(modified_base), rev_strand(is_rev) {}
  mod_id_code() : canonical(' '), value(0), rev_strand(0) {}

  // Get this modification as forward strand mod, i.e. A+a as A+a and T-a as A+a
  mod_id_code const get_plus() const {
    if (this->rev_strand) {
      return mod_id_code(complement(canonical), value, false);
    } else {
      return *this;
    }
  }
  const bool is_plus() { return !this->rev_strand; }

  // Return a copy of this with undefined modification on this kind of
  // nucleotide
  mod_id_code get_undef() { return mod_id_code(canonical, 0, 0); }

  // Check if the modification code is undefined.
  bool is_undef() const { return (value == 0); }

  // Get string representation
  std::string str() const {
    std::stringstream s;
    s << canonical << (rev_strand ? '-' : '+');
    ;
    if (value <= 0) {
      s << -value;
    } else {
      s << std::string(1, (char)value);
    }
    return s.str();
  }



  // Needed for comparisons
  bool operator<(const mod_id_code &other) const {
    return value < other.value ||
           (value == other.value && canonical < other.canonical) ||
           (value == other.value && canonical == other.canonical &&
            rev_strand < other.rev_strand);
  }
  bool operator==(const mod_id_code &other) const {
    return value == other.value && canonical == other.canonical &&
           rev_strand == other.rev_strand;
  }
};

// All mods found in the scanned area
std::set<mod_id_code> found_mods;


class Modification {
public:
  std::string def_str; // Full definition string for this modification type,
                       // context, strand, modification, base position e.g.
                       // "GC+m.1"  or 'T-a.0'
  std::string mod_str; // Character symbol, or numeric ID for this chemical
                       // modifictation (Ignoring base context) e.g. 'm' or 'a'
  std::string fwd_context; // Sequence context in reference forward strand for
                           // reads mapping in forward strand, e.g. "GC" or 'A'
  std::string rev_context; // Sequence context in reference forward strand for
                           // reads mapping in reverse strand, e.g. "GC" or 'T'
  char canonical; // Read original strand base for the modification, e.g. 'C' or
                  // 'T'

  mod_id_code mod_code; // Numerical code matching mod_str
  int fwd_ctx_pos; // Position of the modified base in the fwd_context, e.g. 1
                   // or 0
  int rev_ctx_pos;
  int rev_strand; // Strand of the modification with respect to the original
                  // read strand, e.g. 0 for '+' strand and 1 for '-' strand
  int missing_is_unmodified;

  // E.g.  For methylation modification on C:s at CG context:
  // Modification("CG+m.0")
  Modification(const char *def_str) {
    this->def_str = std::string(def_str);
    // this->mod_code = (mod_id_code)std::hash<std::string>{}(this->def_str);
    std::cmatch cm;
    std::regex mod_pat("([ACGTN]+)([+-])([0-9a-z])([.?])([0-9+])");

    if (!std::regex_match(def_str, cm, mod_pat)) {
      std::cerr << "Couldn't understand modification code '" << def_str
                << std::endl;
      exit(1);
    };

    rev_strand = (cm[2].str()[0] == '-');
    // if (rev_strand) {
    //   std::cerr << "Sorry, can't handle reverse strand modifications"
    //             << std::endl;
    //   // exit(1);
    // }

    mod_str = std::string(cm[3].str());
    char *p;
    int modified_base = -strtol(mod_str.c_str(), &p, 10);
    if (*p) {
      if (mod_str.length() == 1) // require single base modification code
      {
        modified_base = mod_str.c_str()[0];
      } else {
        std::cerr << "Couldn't parse modification type '" << mod_str << "'\n";
        exit(1);
      }
    }
    fwd_context = std::string(cm[1]);
    rev_context = revComplement(fwd_context);

    // fwd_ctx_pos = atoi(cm[5].str().c_str());
    fwd_ctx_pos = strtol(cm[5].str().c_str(), &p, 10);
    rev_ctx_pos = fwd_context.length() - fwd_ctx_pos - 1;
    if (rev_strand) {

      auto tmp = fwd_ctx_pos;
      fwd_ctx_pos = rev_ctx_pos;
      rev_ctx_pos = tmp;
    }

    assert(fwd_ctx_pos >= 0);
    assert(rev_ctx_pos >= 0);

    canonical = fwd_context[fwd_ctx_pos];

    assert(canonical == complement(rev_context[rev_ctx_pos]));

    this->mod_code =
        mod_id_code(canonical, modified_base, this->rev_strand != 0);

    missing_is_unmodified = (cm[4].str()[0] == '.');

    if (!missing_is_unmodified) {
      std::cerr << "Sorry, can't handle ambiguous modification calls"
                << std::endl;
      exit(1);
    }
    // rev_context ==
    //:([ACGTUN][-+]([a-z]+|[0-9]+)[.?]?
    // std::string s = std::string(def);
    // s.find("+")
  }
  std::string to_string() {
    std::stringstream s;
    s << fwd_context << " " << mod_str << canonical << " " << fwd_ctx_pos << ' '
      << (rev_strand ? '-' : '+') << '\n'
      << rev_context << " " << mod_str << canonical << " " << rev_ctx_pos << ' '
      << (rev_strand ? '+' : '-');
    return s.str();
  };
  // True iff modification motif is palindromic.
  bool is_palindromic() { return (fwd_context == rev_context); }
};

std::vector<Modification> modifications;

// Generic stream operator for any map
template <typename V>
std::ostream &operator<<(std::ostream &os, const std::map<mod_id_code, V> &m) {
  os << "{";
  bool first = true;
  for (typename std::map<mod_id_code, V>::const_iterator it = m.begin();
       it != m.end(); ++it) {
    if (!first)
      os << ", ";
    os << it->first.str() << ": " << it->second;
    first = false;
  }
  os << "}";
  return os;
}

class _mod_count_t {
public:
  int total_count; // Total, excluding "other", i.e. canonical + uncalled
                   // sum(modified_count) == total

  int other_count;

  std::map<mod_id_code, int> modified_count;
  std::map<mod_id_code, int> uncalled_count;
  std::map<mod_id_code, double> modified_expected;
  // Overload the stream insertion operator
  friend std::ostream &operator<<(std::ostream &os, const _mod_count_t &p) {
    return os << "total_count:" << p.total_count
              << ", other_count: " << p.other_count
              << "\nmodified_count:" << p.modified_count
              << "\nuncalled_count:" << p.uncalled_count
              << "\nmodified_expected:" << p.modified_expected << '\n';
  }

  _mod_count_t() {
    other_count = 0;
    total_count = 0;
    for (std::vector<Modification>::iterator cmod = modifications.begin();
         cmod != modifications.end(); cmod++) {

      const mod_id_code count_code =
          (combine_hemi ? cmod->mod_code.get_plus() : cmod->mod_code);

      modified_count[count_code] = 0;
      uncalled_count[count_code] = 0;
      modified_expected[count_code] = 0.;
    }
  }
  bool operator==(const _mod_count_t &Ref) const {
    return (this->uncalled_count == Ref.uncalled_count) &&
           this->other_count == Ref.other_count &&
           this->modified_count == Ref.modified_count;
  }
  _mod_count_t &operator+=(_mod_count_t &c) {

    this->other_count += c.other_count;
    this->total_count += c.total_count;

    std::set<mod_id_code> reported_mods{};

    for (std::vector<Modification>::iterator cmod = modifications.begin();
         cmod != modifications.end(); cmod++) {
      const mod_id_code count_code =
          (combine_hemi ? cmod->mod_code.get_plus() : cmod->mod_code);
      if (reported_mods.find(count_code) ==
          reported_mods
              .end()) { // Only sum matching modifications, not contexts
        reported_mods.insert(count_code);

        this->modified_count[count_code] += c.modified_count.at(count_code);
        this->modified_expected[count_code] +=
            c.modified_expected.at(count_code);
        this->uncalled_count[count_code] += c.uncalled_count.at(count_code);
      }

      // std::cerr << this->modified_count.at(cmod->mod_code) << '+' <<
      // c.modified_count.at(cmod->mod_code) << '=' <<
      // this->modified_count[cmod->mod_code] << '@' << cmod->mod_code << '\n';
    }
    return *this;
  }
  _mod_count_t operator+(_mod_count_t &c) {
    _mod_count_t r;

    r.other_count = c.other_count + this->other_count;
    r.total_count = c.total_count + this->total_count;
    for (std::vector<Modification>::iterator cmod = modifications.begin();
         cmod != modifications.end(); cmod++) {
      const mod_id_code count_code =
          (combine_hemi ? cmod->mod_code.get_plus() : cmod->mod_code);
      r.modified_count[count_code] =
          this->modified_count.at(count_code) + c.modified_count.at(count_code);
      r.modified_expected[count_code] = this->modified_expected.at(count_code) +
                                        c.modified_expected.at(count_code);
      r.uncalled_count[count_code] =
          this->uncalled_count.at(count_code) + c.uncalled_count.at(count_code);

      // std::cerr << this->modified_count.at(cmod->mod_code) << '+' <<
      // c.modified_count.at(cmod->mod_code) << '=' <<
      // r.modified_count[cmod->mod_code] << '@' << cmod->mod_code << '\n';
    }
    return r;
  }

  void add_uncalled(mod_id_code mod_id) { uncalled_count[mod_id]++; }
  void add_match() { total_count++; }
  void add_modified(mod_id_code mod_id) {
    assert((!combine_hemi) || mod_id.is_plus());
    modified_count[mod_id]++;
    modified_expected[mod_id] += 1.;
  }
  /* Add expectation and count modification for give mod.
   */
  void add_observed(mod_id_code mod_id, double mod_expectation,
                    bool is_modified) {

    if (is_modified) {
      modified_count[mod_id]++;
    }
    modified_expected[mod_id] += mod_expectation;
  }
  void add_other() { other_count++; }
  int get_called_count(const mod_id_code &mod_id) {
    return this->total_count - this->uncalled_count.at(mod_id);
  }
  int get_uncalled_count(const mod_id_code &mod_id) {
    return this->uncalled_count.at(mod_id);
  }
  int get_modified_count(mod_id_code mod_id) {
    return this->modified_count.at(mod_id);
  }
};

std::ostream &operator<<(std::ostream &ss, _mod_count_t &obj) {
  ss << obj.total_count << '\t' << obj.other_count;
  for (std::map<mod_id_code, int>::iterator it = obj.modified_count.begin();
       it != obj.modified_count.end(); it++) {
    ss << '\t';
    ss << it->first.str() << ":" << it->second << ":"
       << obj.uncalled_count.at(it->first);
  }
  return ss;
}

class mod_key {
public:
  mod_key(int pos = 0, int read_is_rev = -1, int ps = -1, int hp = -1)
      : pos(pos), read_isrev(read_is_rev), ps(ps), hp(hp) {}
  bool operator==(const mod_key &Ref) const {
    return (this->pos == Ref.pos) && (this->read_isrev == Ref.read_isrev) &&
           (this->ps == Ref.ps) && (this->hp == Ref.hp);
  }

  // Order of mod_keys: dictionary order of (pos,read_is_rev,ps,hp)
  bool operator<(const mod_key &rhs) const {
    if (pos < rhs.pos)
      return true;
    else if (pos > rhs.pos)
      return false;
    if (read_isrev < rhs.read_isrev)
      return true;
    else if (read_isrev > rhs.read_isrev)
      return false;
    if (ps < rhs.ps)
      return true;
    else if (ps > rhs.ps)
      return false;
    if (hp < rhs.hp)
      return true;
    else if (hp > rhs.hp)
      return false;

    return false;
  }
  const int pos;

  const int read_isrev;
  const int ps;
  const int hp;
};
std::ostream &operator<<(std::ostream &ss, const mod_key &obj_) {
  ss << obj_.pos << '\t' << obj_.read_isrev << '\t' << obj_.ps << '\t'
     << obj_.hp;
  return ss;
}

class ModCounter {
private:
  /* data */
  std::map<mod_key, _mod_count_t> mod_count_store;

public:
  ModCounter(/* args */);
  ~ModCounter();
  void add_uncalled(int pos, int read_is_rev, int ps, int hp,
                    mod_id_code mod_id);
  void add_canonical_except(int pos, int read_is_rev, int ps, int hp,
                            std::set<mod_id_code> mods_id_rev);
  void add_modified(int pos, int read_is_rev, int ps, int hp,
                    mod_id_code mod_id, double mod_expect = 1.,
                    bool is_modified = true);

  void add_other(int pos, int read_is_rev, int ps = -1, int hp = -1);
  void add_match(int pos, int read_is_rev, int ps = -1, int hp = -1) {
    mod_key k(pos, read_is_rev, ps, hp);
    mod_count_store[k].add_match();
  }
  // For debugging only
  _mod_count_t return_inner(int pos, int read_is_rev, int ps = -1,
                            int hp = -1) {
    mod_key k(pos, read_is_rev, ps, hp);
    return mod_count_store.at(k);
  }
  std::map<mod_key, _mod_count_t>::iterator begin() {
    return mod_count_store.begin();
  }
  std::map<mod_key, _mod_count_t>::iterator end() {
    return mod_count_store.end();
  }
  std::size_t size() { return this->mod_count_store.size(); }
  std::map<mod_key, _mod_count_t>::iterator from(int pos, int read_is_rev = -1,
                                                 int ps = -1, int hp = -1) {
    std::map<mod_key, _mod_count_t>::iterator it =
        mod_count_store.lower_bound(mod_key(pos, read_is_rev, ps, hp));
    // if (it != this->end())
    // {
    //     std::cout << this->begin()->second << ' ';
    //     std::cout << pos << '=' << it->first << "+" << (*it).second << '\n';
    // }

    // assert(it != this->end());
    return it;
  }
  std::string output_mod_joinstrand(char const *ref_name, int pos,
                                    Modification &cmod,
                                    WhichStrandOutput which_strand);
  void check() {
    if (!mod_count_store.empty()) {
      // std::cerr << mod_count_store.begin()->first << '@' <<
      // mod_count_store.begin()->second << " <-> " <<
      // mod_count_store.at(mod_count_store.begin()->first) << '\n';
      assert(mod_count_store.count(mod_count_store.begin()->first) == 1);
      assert(mod_count_store.begin()->second ==
             mod_count_store.at(mod_count_store.begin()->first));
    }
    for (std::map<mod_key, _mod_count_t>::iterator it = mod_count_store.begin();
         it != mod_count_store.end(); it++) {
      // std::cerr << mod_count_store.at(it->first) << " --- " << it->second
      //           << " ... " << mod_count_store[it->first] << '\n';
      assert(mod_count_store.at(it->first) == it->second);
      assert(mod_count_store.find(it->first)->first == it->first);
#ifndef NDEBUG
      std::map<mod_key, _mod_count_t>::iterator fit =
          mod_count_store.lower_bound(it->first);
#endif
      assert(fit == it);
      assert(fit->first == it->first);
      assert(fit->second == it->second);
    }
    // if (mod_count_store.size() > 0)
    // abort();
  }
  std::map<mod_key, _mod_count_t>::iterator to(int pos) {
    return mod_count_store.lower_bound(mod_key(pos + 1));
  }
  // std::map<mod_key, _mod_count_t>::iterator at(int pos) { return
  // mod_count_store.find(mod_key(pos)); }
  void clear() { mod_count_store.clear(); }
  void erase_upto(int upto);
};

std::ostream &operator<<(std::ostream &os, ModCounter &cnts) {
  for (std::map<mod_key, _mod_count_t>::iterator it = cnts.begin();
       it != cnts.end(); it++) {
    os << it->first << " : " << it->second << '\n';
  }
  return os;
}
ModCounter::ModCounter(/* args */) {
  // this->mod_count_store = std::map<mod_key, _mod_count_t>();
}

ModCounter::~ModCounter() {}

ModCounter mod_counter;

// Remove all modifications up to, and including position 'upto'
void ModCounter::erase_upto(int upto) {
  if (mod_count_store.size() > 0) {

    std::map<mod_key, _mod_count_t>::iterator upto_it =
        this->mod_count_store.lower_bound(upto + 1);

    this->mod_count_store.erase(this->begin(), upto_it);
  }
}

std::string mod_count_to_str(mod_id_code mod_code, int mod_count,
                             int called_sites) {
  std::stringstream ss;
  ss << mod_code.str() << '\t' << mod_count;
  double mod_freq = (double)mod_count / (double)called_sites;
  ss << '\t' << mod_freq;
  return ss.str();
}

// Motif match strand: which_strand == 0 forward, == 1 reverse, -1==both
std::string ModCounter::output_mod_joinstrand(char const *chrom, int pos,
                                              Modification &cmod,
                                              WhichStrandOutput which_strand) {
  std::stringstream ss, out_stream;
  _mod_count_t mod_counts_total;

  if (header_flag) {
    out_stream
        << "#chromosome\tstart\tend\thaplotype\tphase_set\tcalled_"
           "reads\tuncalled_reads\tmismatch_"
           "reads\tmodification\tmodified_reads\tmodified_prop\tcall_code\n";
  }

  ss << chrom << '\t' << pos + 2 - cmod.fwd_context.length() << '\t' << pos + 2
     << '\t';

  int fwd_pos = pos - cmod.fwd_context.length() + cmod.fwd_ctx_pos + 1;
  int rev_pos = pos - cmod.fwd_context.length() + cmod.rev_ctx_pos + 1;

  // pair<phase_set,haplotype> -> _mod_count_t
  std::map<std::pair<int, int>, _mod_count_t> cnts;

  if (fwd_pos == rev_pos) {
    // Sum the modifications over both read alignment strands

    for (std::map<mod_key, _mod_count_t>::iterator it = this->from(fwd_pos);
         it != this->to(fwd_pos); it++) {

      cnts[std::pair<int, int>(it->first.ps, it->first.hp)] =
          cnts[std::pair<int, int>(it->first.ps, it->first.hp)] + it->second;
    }
  } else {

    // If the motif has an other modification location in one of the reference
    // strands, sum those in correct strands
    int first_pos = (fwd_pos < rev_pos ? fwd_pos : rev_pos);
    int last_pos = (fwd_pos > rev_pos ? fwd_pos : rev_pos);

    for (std::map<mod_key, _mod_count_t>::iterator it = this->from(first_pos);
         it != this->to(last_pos); it++) {
      if ((it->first.read_isrev == 0 && it->first.pos == fwd_pos) ||
          (it->first.read_isrev == 1 && it->first.pos == rev_pos)) {
        cnts[std::pair<int, int>(it->first.ps, it->first.hp)] =
            cnts[std::pair<int, int>(it->first.ps, it->first.hp)] + it->second;
      }
    }
  }

  // if (which_strand == STRAND_FORWARD || which_strand == STRAND_BOTH) {
  //   // std::cerr << "pos : " << pos << " +" << fwd_pos << " -" << rev_pos <<
  //   // '\n';
  //   for (std::map<mod_key, _mod_count_t>::iterator it = this->from(fwd_pos,
  //   0);
  //        it != this->to(fwd_pos); it++) {
  //     // std::cerr << "fwd it: " << it->first << " -  " << it->second <<
  //     '\n'; if (it->first.read_isrev == 0) {
  //       // std::cerr << it->first.ps << ',' << it->first.hp << " +@ " <<
  //       // cnts[std::pair(it->first.ps, it->first.hp)] << " + " <<
  //       it->second;

  //       cnts[std::pair<int, int>(it->first.ps, it->first.hp)] =
  //           cnts[std::pair<int, int>(it->first.ps, it->first.hp)] +
  //           it->second;
  //       // std::cerr << " = " << cnts[std::pair(it->first.ps, it->first.hp)]
  //       <<
  //       // '\n';
  //     }
  //   }
  // }

  // if (which_strand == STRAND_REVERSE || which_strand == STRAND_BOTH) {
  //   for (std::map<mod_key, _mod_count_t>::iterator it = this->from(rev_pos,
  //   0);
  //        it != this->to(rev_pos); it++) {
  //     if (it->first.read_isrev == 1) {
  //       // std::cerr << it->first.ps << ',' << it->first.hp << " -@ " <<
  //       // cnts[std::pair(it->first.ps, it->first.hp)] << " + " <<
  //       it->second; cnts[std::pair<int, int>(it->first.ps, it->first.hp)] =
  //           cnts[std::pair<int, int>(it->first.ps, it->first.hp)] +
  //           it->second;
  //       // std::cerr << " = " << cnts[std::pair(it->first.ps, it->first.hp)]
  //       <<
  //       // '\n';
  //     }
  //   }
  // }

  // Ensure combining the A+a and T-a when needed.
  const mod_id_code mod_code = cmod.mod_code;
  const mod_id_code counted_mod_code =
      (combine_hemi ? mod_code.get_plus() : mod_code);
  //(combine_hemi ? cmod.mod_code.get_plus() : cmod.mod_code);

  // Print phase set and haplotype totals:
  for (std::map<std::pair<int, int>, _mod_count_t>::iterator cnts_it =
           cnts.begin();
       cnts_it != cnts.end(); cnts_it++) {
    const int haplotype = cnts_it->first.second;
    const int phase_set = cnts_it->first.first;

    _mod_count_t &mod_counts_here = cnts_it->second;

    const int called_sites = mod_counts_here.get_called_count(counted_mod_code);
    const int modified_sites =
        mod_counts_here.get_modified_count(counted_mod_code);
    const int uncalled_sites =
        mod_counts_here.get_uncalled_count(counted_mod_code);

    assert(called_sites <= mod_counts_here.total_count);
    out_stream << ss.str();
    if (haplotype >= 0) {
      out_stream << haplotype << '\t' << phase_set << '\t';
    } else {
      out_stream << "N\t-1\t";
    }
    assert(mod_counts_here.total_count >= 0);
    assert(called_sites >= 0);
    assert(mod_counts_here.total_count >= called_sites);
    assert(mod_counts_here.other_count >= 0);
    assert(mod_counts_here.modified_count[mod_code] >= 0);

    out_stream << called_sites << '\t' << uncalled_sites << '\t'
               << mod_counts_here.other_count << '\t';

    out_stream << mod_count_to_str(mod_code, modified_sites, called_sites)
               << '\t' << cmod.def_str << '\n';

    mod_counts_total += mod_counts_here;
  }
  // Print the locus total:
  const int called_sites = mod_counts_total.get_called_count(counted_mod_code);
  const int modified_sites =
      mod_counts_total.get_modified_count(counted_mod_code);
  const int uncalled_sites =
      mod_counts_total.get_uncalled_count(counted_mod_code);
  assert(called_sites <= mod_counts_total.total_count);
  assert(mod_counts_total.total_count >= 0);
  assert(called_sites >= 0);
  assert(mod_counts_total.total_count >= called_sites);
  assert(mod_counts_total.other_count >= 0);
  assert(mod_counts_total.modified_count[mod_code] >= 0);

  if (called_sites == 0) {
    return std::string("");
  } else {
    out_stream << ss.str() << "*\t*\t" << called_sites << '\t' << uncalled_sites
               << '\t' << mod_counts_total.other_count << '\t';
    out_stream << mod_count_to_str(mod_code, modified_sites, called_sites)
               << '\t' << cmod.def_str;

    out_stream << '\n';
    header_flag = 0;
    return out_stream.str();
  }
}
void ModCounter::add_uncalled(int pos, int read_is_rev, int ps, int hp,
                              mod_id_code mod_id) {
  mod_key k(pos, read_is_rev, ps, hp);

  auto &mod_cnts = mod_count_store[k];

  mod_cnts.add_uncalled(mod_id);
}

void ModCounter::add_other(int pos, int read_is_rev, int ps, int hp) {
  mod_key k(pos, read_is_rev, ps, hp);
  mod_count_store[k].add_other();
}

void ModCounter::add_modified(int pos, int read_is_rev, int ps, int hp,
                              mod_id_code mod_id, double mod_expect,
                              bool is_modified) {
  // assert(mod_is_rev == 0);
  mod_key k(pos, read_is_rev, ps, hp);
  mod_count_store[k].add_observed(mod_id, mod_expect, is_modified);
}

void output_mod_counter(int &pos) {

  for (int read_is_rev : {0, 1}) {
    for (int hp : {1, 2}) {
      int PS = 45553715;
      std::cerr << '\n'
                << pos << ":isRev" << read_is_rev << ":PS" << PS << ":haplo"
                << hp << '\n'
                << mod_counter.return_inner(pos, read_is_rev, PS, hp);
    }
  }
}

#define MAX_MOD_COUNT 10
/* Report a line of pileup, including base modifications inline with
 the sequence (including qualities), as [<strand><dir><qual>...]

 Counting here the modifications per nucleotide, strand and modification (like
 'm' or 'a'). This does not involve with sequence context or anything else.

 Filtering is by min_baseq.
 */

void process_mod_pileup0(sam_hdr_t *h, const bam_pileup1_t *p, int tid, int pos,
                         int n, char ref_base) {

  int i;

  for (i = 0; i < n; i++, p++) {
    uint8_t *seq = bam_get_seq(p->b);
    uint8_t *qual = bam_get_qual(p->b);
    unsigned char q_base = seq_nt16_str[bam_seqi(seq, p->qpos)];
    uint8_t base_qual = qual[p->qpos];
    int ps = -1, hp = -1;
    if (base_qual < min_baseq) { // Ignore bases with bad quality.

      continue;
    }
    int read_is_rev = bam_is_rev(p->b);

    if (phase_flag) {
      uint8_t *_hp_tag = bam_aux_get(p->b, "HP");

      if (_hp_tag != NULL) {
        hp = bam_aux2i(_hp_tag);
        ps = bam_aux2i(bam_aux_get(p->b, "PS"));
      }
    }

    if (p->is_del || q_base != ref_base) {
      mod_counter.add_other(pos, read_is_rev, ps,
                            hp); // Read vs reference base mismatch
    } else {
      mod_counter.add_match(pos, read_is_rev, ps, hp);
      hts_base_mod_state *m = (hts_base_mod_state *)p->cd.p;

      hts_base_mod mod[MAX_MOD_COUNT];
      int nm;
      if ((nm = bam_mods_at_qpos(p->b, p->qpos, m, mod, MAX_MOD_COUNT)) > 0) {
        int j;
        
        for (j = 0; j < nm && j < MAX_MOD_COUNT; j++) {

          assert(
              (!read_is_rev && mod[j].canonical_base == ref_base) ||
              (read_is_rev && mod[j].canonical_base == complement(ref_base)));
          mod_id_code mod_id(mod[j].canonical_base, mod[j].modified_base,
                             mod[j].strand);

          found_mods.insert(mod_id);
          if (combine_hemi) { // If combine hemimethylated calls, i.e. count
                              // A+a and T-a as one
            mod_id = mod_id.get_plus();
          }
          assert(!mod_id.is_undef());
          assert((!combine_hemi) || mod_id.is_plus());

          if ((mod[j].qual > max_mod_prob)) {
            if (mod[j].qual < min_mod_prob) {
              // Modification has ambiguous likelihood
              mod_counter.add_uncalled(pos, read_is_rev, ps, hp, mod_id);
            } else {

              double modified_expected = mod[j].qual / 255.;
              mod_counter.add_modified(pos, read_is_rev, ps, hp, mod_id,
                                       modified_expected);
            }
          }
        }
      }
    }
  }
  // if (pos == 46138257)
  //   output_mod_counter(pos);
}

std::ostream &output_mod(std::ostream &out_stream, const mod_key &mod_id,
                         _mod_count_t &cnts, const char *chrom,
                         Modification &out_mod) {
  // std::cout << mod_id << 'x' << cnts << '\n';

  const mod_id_code count_code =
      (combine_hemi ? out_mod.mod_code.get_plus() : out_mod.mod_code);

  const int called_sites = cnts.get_called_count(count_code);
  const int modified_sites = cnts.get_modified_count(count_code);

  assert(called_sites <= cnts.total_count);
  std::stringstream ss;
  if (header_flag) {
    out_stream
        << "#chromosome\tposition\tstrand\thaplotype\tphase_set\tcalled_"
           "reads\tuncalled_reads\tmismatch_"
           "reads\tmodification\tmodified_reads\tmodified_prop\tcall_code\n";
    header_flag = 0;
  }
  // abort();
  //  for (std::map<int, int>::iterator it = cnts.modified_count.begin();
  //       it != cnts.modified_count.end(); it++)

  // {
  //     called_sites += it->second;
  // }

  ss << chrom << '\t' << mod_id.pos + 1 << '\t'
     << (mod_id.read_isrev ? '-' : '+') << '\t';

  if (mod_id.hp >= 0) {
    ss << mod_id.hp << '\t' << mod_id.ps << '\t';
  } else {
    ss << "N\t" << mod_id.ps << '\t';
  }
  const int uncalled_sites = cnts.get_uncalled_count(count_code);
  ss << called_sites << '\t' << uncalled_sites << '\t' << cnts.other_count
     << '\t';

  out_stream << ss.str();

  out_stream << mod_count_to_str(out_mod.mod_code, modified_sites, called_sites)
             << '\t' << out_mod.def_str;

  return out_stream << '\n';
}

static char *reference_fasta_file = NULL;
static char *input_bam_file = NULL;
static char *target_region = NULL;
#define PARSE_UNKNOWN 1
#define PARSE_HELP 2
#define PARSE_VERSION 3

int parse_options(int argc, char **argv) {
  int c;
  static struct option long_options[] = {
      /* These options set a flag. */
      {"no_header", no_argument, &header_flag, 0},
      {"no_phase", no_argument, &phase_flag, 0},
      {"split_strand", no_argument, &split_strand_flag, 1},
      {"no_combine_mod_strands", no_argument, (int *)&combine_hemi, 0},
      /* These options don’t set a flag.
         We distinguish them by their indices. */
      {"reference_fasta", required_argument, 0, 'r'},
      {"input", required_argument, 0, 'i'},
      {"min_mod_prob", required_argument, 0, 'c'},
      {"max_mod_prob", required_argument, 0, 'C'},
      {"min_baseq", required_argument, 0, 'b'},
      {"min_mapq", required_argument, 0, 'q'},
      {"exclude", required_argument, 0, 'E'},
      {"mod", required_argument, 0, 'm'},
      {"region", required_argument, 0, 'R'},
      {"version", no_argument, 0, 'V'},
      {"help", no_argument, 0, 'h'},
      {0, 0, 0, 0}};
  while (1) {

    /* getopt_long stores the option index here. */
    int option_index = 0;

    c = getopt_long(argc, argv, "r:i:c:C:b:m:q:E:R:@:hV", long_options,
                    &option_index);
    // std::cerr << "getopt_long ret: " << c << '\n';
    /* Detect the end of the options. */
    if (c == -1)
      break;

    // std::cerr << "phase_flag:" << phase_flag << '\n';
    if (c != 0) {
      std::cerr << "option -" << (char)c << ' ' << (optarg ? optarg : "")
                << '\n';
    } // std::cerr << c << ':' << option_index << ':' <<
      // long_options[option_index].name << ':' << optarg << std::endl;
    // std::cerr << "Err" << '\n';
    // std::cout << __LINE__ << '\n';
    switch (c) {
    case 0:
      /* If this option set a flag, do nothing else now. */
      if (long_options[option_index].flag != 0)
        break;
      std::cerr << "option --" << long_options[option_index].name;
      if (optarg)
        std::cerr << optarg;
      std::cerr << '\n';
      break;

    case 'm':

      modifications.push_back(Modification(optarg));
      break;
    case 'r':
      // printf("option -r with value `%s'\n", optarg);
      reference_fasta_file = strdup(optarg);
      break;

    case 'i':
      // printf("option -i with value `%s'\n", optarg);
      input_bam_file = strdup(optarg);
      break;

    case 'c':
      // printf("option -c with value `%s'\n", optarg);
      min_mod_prob = atoi(optarg);
      break;
    case 'C':
      max_mod_prob = atoi(optarg);
      break;
    case 'b':
      // printf("option -b with value `%s'\n", optarg);
      min_baseq = atoi(optarg);
      break;
    case 'q':
      // printf("option -q with value `%s'\n", optarg);
      min_mapq = atoi(optarg);
      break;
    case 'E':
      // printf("option -E with value `%s'\n", optarg);
      exclude_filter = atoi(optarg);
      break;
    case 'R':
      // printf("option -R with value `%s'\n", optarg);
      target_region = strdup(optarg);
      break;
    case '@':
      extra_hts_threads = atoi(optarg);
      break;
    case 'h':
      return PARSE_HELP;
    case 'V':
      fprintf(stderr, "Version %d.%d.%d\n", btm_version.major,
              btm_version.minor, btm_version.patch);
      return PARSE_VERSION;
    case '?':
      /* getopt_long already printed an error message. */
      return PARSE_UNKNOWN;
      break;

    default:
      abort();
    }
  }
  if (modifications.size() == 0) {
    modifications.push_back(Modification("CG+m.0"));
  }

  if (!reference_fasta_file || !input_bam_file) {
    return 1;
  }
  return 0;
}
/**
 * @brief Change all characters in the string arr to upper case.
 *
 * @param arr
 * @return char*
 */
char *to_upper_char_arr(char *arr) {
  char *p;
  for (p = arr; *p; p++) {
    *p = (char)toupper(*p);
  }
  return arr;
}
#if HTS_VERSION < 101600
#error "Requires HTSLIB version at least 1.16. "
#endif
int main(int argc, char **argv) {
  if (int r = parse_options(argc, argv)) {
    if (r == PARSE_HELP) {

      fprintf(stderr,
              "usage: %s [-c %d] [-C %d]  [-b %d]  [-q %d] [-E %d] [-m CG+m.0] "
              "[-R chr1:1-100] [--no_header] [--no_phase] [--split_strand] "
              "[--no_combine_mod_strands] -r "
              "ref.fasta -i input.cram\n\n"
              " -c|--min_mod_prob Minimum probability of modification called "
              "modified (range 0-255)\n"
              " -C|--max_mod_prob Maximum probability of modification called "
              "unmodified (range 0-255)\n"
              " -b|--min_baseq    Minimum base quality considered.\n"
              " -q|--min_mapq     Minimum mapping quality considered.\n"
              " -E|--exclude      Exclude all reads matching any of these SAM "
              "flags.\n"
              " -R|--region       Genomic region to consider.\n"
              " -r|--reference_fasta   FAI indexed fasta file of the reference "
              "genome.\n"
              " -i|--input        Input SAM/BAM/CRAM file. Indexed if used "
              "with -R.\n"
              " -m|--mod          DNA modification to consider. Format: "
              "'CG+m.0' for methylation 'm' of C:s on position '0' \n"
              "                   of context 'CG' in forward strand. Can be "
              "given multiple times. If none given, 'CG+m.0' is used.\n"
              " --split_strand    Report modifications on either strand "
              "separately.\n"
              " --no_combine_mod_strands  Don't count hemimethylated loci. "
              "i.e. Count e.g. A+a and T-a as sparate.\n"
              " --no_header       Don't output header text.\n"
              " --no_phase        Don't use HP and PS tags to separate phased "
              "reads.\n"
              " -@                Extra threads for reading input.\n"
              " -h|--help         Print this help text.\n\n"
              "Version %d.%d.%d\n"
              "Using htslib version %s.\n",
              argv[0], min_mod_prob, max_mod_prob, min_baseq, min_mapq,
              exclude_filter, btm_version.major, btm_version.minor,
              btm_version.patch, hts_version());
    }

    // Don't fail if asked for help.

    exit((r == PARSE_HELP || r == PARSE_VERSION) ? 0 : 1);
  }

  if (strcmp(hts_version(), "1.16") < 0) {
    fprintf(stderr,
            "Need htslib version at least 1.16. Compiled with %d. Currently "
            "having %s.\n",
            HTS_VERSION, hts_version());
    exit(2);
  }
  if (strcmp(hts_version(), "1.16") < 0) {
    fprintf(stderr,
            "Need htslib version at least 1.16. Compiled with %d. Currently "
            "having %s.\n",
            HTS_VERSION, hts_version());
    exit(2);
  }

  // if(split_strand_flag) {
  //   std::cerr<<"Reporting modifications by read alignment strand. Reporting
  //   only exact strand mods, i.e. A+a and T-a are not merged.\n";
  //   combine_hemi=false;
  // }

  // First argument: reference genome
  std::cerr << "reference_fasta_file:" << reference_fasta_file << std::endl;
  std::cerr << "input_bam_file:" << input_bam_file << std::endl;
  std::cerr << "min_mod_prob:" << min_mod_prob << " = "
            << (min_mod_prob + 0.5) / 256. << std::endl;
  std::cerr << "max_mod_prob:" << max_mod_prob << " = "
            << (max_mod_prob + 0.5) / 256. << std::endl;
  std::cerr << "Version: " << btm_version.major << '.' << btm_version.minor
            << '.' << btm_version.patch << std::endl;
  // First argument: reference genome
  std::cerr << "reference_fasta_file:" << reference_fasta_file << std::endl;
  std::cerr << "input_bam_file:" << input_bam_file << std::endl;
  std::cerr << "min_mod_prob:" << min_mod_prob << " = "
            << (min_mod_prob + 0.5) / 256. << std::endl;
  std::cerr << "max_mod_prob:" << max_mod_prob << " = "
            << (max_mod_prob + 0.5) / 256. << std::endl;
  std::cerr << "split_strands:" << (split_strand_flag ? "True\n" : "False\n");
  std::cerr << "Version: " << btm_version.major << '.' << btm_version.minor
            << '.' << btm_version.patch << std::endl;

  for (std::vector<Modification>::iterator itr = modifications.begin();
       itr != modifications.end(); itr++) {
    std::cerr << itr->to_string() << std::endl;
  }

  faidx_t *ref_fai = fai_load(reference_fasta_file);
  if (!ref_fai) {
    fprintf(stderr, "Can't open reference genome fasta '%s'.",
            reference_fasta_file);

    exit(1);
  }

  samFile *in = sam_open(input_bam_file, "r");
  if (!in) {
    fprintf(stderr, "Can't open input bam file '%s'.", input_bam_file);

    exit(1);
  }

  if (extra_hts_threads > 0) {
    hts_set_threads(in, extra_hts_threads);
  }

  bam1_t *b = bam_init1();
  sam_hdr_t *h = sam_hdr_read(in);

  const bam_pileup1_t *p;
  int tid, pos, n;
  hts_pos_t begin = 0, end = INT_MAX;
  int prev_tid = -1;
  hts_pos_t ref_len = 0;
  char *ref_seq = NULL;
  const char *ref_seq_name = NULL;

  hts_itr_t *sam_iter = NULL;
  if (target_region != NULL) {
    sam_parse_region(h, target_region, &tid, &begin, &end, 0);

    hts_idx_t *idx = NULL;
    if ((idx = sam_index_load(in, input_bam_file)) == 0) {
      fprintf(stderr, "[E::%s] fail to load the BAM index\n", __func__);
      exit(1);
    }
    if ((sam_iter = sam_itr_querys(idx, h, target_region)) == 0) {
      fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__,
              target_region);
      exit(1);
    }
  }
  // Pileup iterator with constructor/destructor to parse base mod tags
  plp_dat dat = {
      .fp = in,
      .fp_hdr = h,
      .itr = sam_iter,
  };
  bam_plp_t iter = bam_plp_init(readaln, &dat);
  bam_plp_constructor(iter, pileup_cd_create);
  bam_plp_destructor(iter, pileup_cd_destroy);

  int motif_sites_covered = 0;

  while ((p = bam_plp_auto(iter, &tid, &pos, &n)) !=
         0) { // For each reference position
    bool plp_procced = false;

    if (pos > end) {
      break;
    }

    // New contig/chromosome
    if (tid != prev_tid) {
      mod_counter.clear();
      if (ref_seq) {
        free(ref_seq);
      }
      ref_seq_name = sam_hdr_tid2name(dat.fp_hdr, tid);
      ref_seq = fai_fetch64(ref_fai, ref_seq_name, &ref_len);
      if (!ref_seq) {
        fprintf(stderr, "Couldn't fetch reference!");
        exit(1);
      } else {
        std::cerr << "Reading chromosome " << ref_seq_name << '\n';
        ref_seq = to_upper_char_arr(ref_seq);
      }
      prev_tid = tid;
    }
    if (pos >= ref_len) {
      fprintf(stderr, "Trying to access reference beyond end!");
      exit(1);
    }

    // Only output sites with correct context (e.g. CpG sites).
    for (std::vector<Modification>::iterator cmod = modifications.begin();
         cmod != modifications.end();
         cmod++) { // For each queried modification at the reference position
      // require one full context length at beginning and end of reference
      if (pos < (int)cmod->fwd_context.length() ||
          pos > ref_len - (int)cmod->fwd_context.length())
        continue;
      int m_fwd = cmod->fwd_context.compare(0, std::string::npos,
                                            ref_seq + pos - cmod->fwd_ctx_pos,
                                            cmod->fwd_context.length());
      int m_rev = cmod->rev_context.compare(0, std::string::npos,
                                            ref_seq + pos - cmod->rev_ctx_pos,
                                            cmod->fwd_context.length());

      // Correct context
      if (m_fwd == 0 || m_rev == 0) {
        if (!plp_procced) { // Count the modified and unmodified bases.
          process_mod_pileup0(h, p, tid, pos, n, ref_seq[pos]);
          motif_sites_covered++;
          plp_procced = true;
        }
      }

      // Output modification at the given reference position
      // Out of region
      if ((pos <= begin) || (pos > end)) {
        continue;
      }
      if (split_strand_flag) {
        for (std::map<mod_key, _mod_count_t>::iterator it =
                 mod_counter.from(pos);
             it != mod_counter.to(pos); it++) {
          assert(it != mod_counter.end());
          if ((m_fwd == 0 && it->first.read_isrev == 0) ||
              ((m_rev == 0 && it->first.read_isrev == 1))) {

            output_mod(std::cout, it->first, it->second, ref_seq_name, *cmod);
          }
        }
      }

      else {
        // Gone past the 'latter' of the modified sites of the motif.
        if ((cmod->fwd_ctx_pos < cmod->rev_ctx_pos && m_rev == 0) ||
            (cmod->fwd_ctx_pos > cmod->rev_ctx_pos && m_fwd == 0)) {
          std::cout << mod_counter.output_mod_joinstrand(ref_seq_name, pos,
                                                         *cmod, STRAND_BOTH);
        } else if (cmod->fwd_ctx_pos == cmod->rev_ctx_pos &&
                   (m_rev == 0 || m_fwd == 0)) {
          std::cout << mod_counter.output_mod_joinstrand(
              ref_seq_name, pos, *cmod,
              // (combine_hemi ? STRAND_BOTH
              //               : (m_fwd == 0 ? STRAND_FORWARD :
              //               STRAND_REVERSE))
              (m_fwd == 0 ? STRAND_FORWARD : STRAND_REVERSE));
        }
      }
    }

    // Remove the cached methyl counts.
    if (mod_counter.begin()->first.pos < (pos - 100)) {
      mod_counter.erase_upto(pos - 10);
    }
  }
  bam_plp_destroy(iter);

  sam_close(in);
  bam_destroy1(b);
  sam_hdr_destroy(h);
  free(ref_fai);
  std::cerr << "Covered " << motif_sites_covered << " sites.\nEncountered modifications:";

  for(const mod_id_code &mod_id : found_mods){
    std::cerr << ' ' << mod_id.str();
  }
  std::cerr<<std::endl;
  // Fail if no motifs covered.
  return (motif_sites_covered > 0 ? 0 : 1);
}
