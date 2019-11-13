// [[Rcpp::plugins("cpp11")]]

#include<Rcpp.h>
#include<string>
#include<unordered_map>
#include<vector>
#include<iostream>
#include<algorithm>
#include<functional>

typedef std::function<std::string(std::string&, int)> KMER_DECORATOR;

const int HASH_CONST = 233;

const int MOD = 1e9 + 33;

int get_window_length(const Rcpp::IntegerVector& d) {
  int res = 0;
  for(const int& item: d) {
    res += item + 1;
  }
  return res;
}

int get_hash(const std::vector<int>& s,
             const Rcpp::IntegerVector& d,
             int begin_index,
             bool pos) {
  long long lres = ((pos ? begin_index : 0) * HASH_CONST + s[begin_index]) % MOD;
  int current_letter_index = begin_index;
  for(int i = 0; i < d.size(); ++i) {
    current_letter_index += d[i] + 1;
    lres = ((lres * HASH_CONST) + s[current_letter_index]) % MOD;
  }
  return (int)lres;
}

int get_hash(const std::vector<int>& kmer) {
  long long hash = 0;
  for(int i = 0; i < kmer.size(); ++i) {
    hash = ((hash * HASH_CONST) + kmer[i]) % MOD;
  }
  return (int)hash;
}

int get_total_size_of_kmer(const std::vector<int>& s,
                           const Rcpp::IntegerVector& d,
                           int begin_index,
                           std::unordered_map<int, std::string>& num2str) {
  int res = num2str[s[begin_index]].size();
  int current_index = begin_index;
  for(int d_i = 0; d_i < d.size(); ++d_i) {
    current_index += d[d_i] + 1;
    res += num2str[s[current_index]].size() + 1; // + 1 because of the dot separator
  }
  return res;
}

int get_total_size_of_kmer(const std::vector<int>& kmer,
                           std::unordered_map<int, std::string>& num2str) {
  int total_size = 0;
  for(const int& item: kmer) {
    total_size += num2str[item].size() + 1;
  }
  return total_size - 1;
}

std::string create_kmer(const std::vector<int>& s,
                        const Rcpp::IntegerVector& d,
                        int begin_index,
                        std::unordered_map<int, std::string>& num2str,
                        KMER_DECORATOR kmer_decorator) {
  int total_kmer_size = get_total_size_of_kmer(s, d, begin_index, num2str);
  std::string res;
  res.reserve(total_kmer_size);
  res += num2str[s[begin_index]];
  int current_index = begin_index;
  for(int d_i = 0; d_i < d.size(); ++d_i) {
    current_index += d[d_i] + 1;
    res += "." + num2str[s[current_index]];
  }
  return kmer_decorator(res, begin_index);
}

std::string create_kmer(const std::vector<int>& kmer,
                        std::unordered_map<int, std::string>& num2str,
                        KMER_DECORATOR kmer_decorator) {
  int total_kmer_size = get_total_size_of_kmer(kmer, num2str);
  std::string res;
  res.reserve(total_kmer_size);
  for(const int& elem: kmer) {
    if(res.size() == 0) {
      res += num2str[elem];
    } else {
      res += "."  + num2str[elem];
    }
  }
  return kmer_decorator(res, -1);
}

void update_kmers(std::unordered_map<int, std::pair<std::string, int>>& kmers,
                  const Rcpp::IntegerVector& d,
                  const std::vector<int>& s,
                  int kmer_hash,
                  int kmer_begin_index,
                  std::unordered_map<int, std::string>& num2str,
                  KMER_DECORATOR kmer_decorator) {
  auto map_elem = kmers.find(kmer_hash);
  if(map_elem == kmers.end()) {
    kmers[kmer_hash] = std::make_pair(create_kmer(s, d, kmer_begin_index, num2str, kmer_decorator), 0);
  }
  ++kmers[kmer_hash].second;
}

void add_kmer_if_not_exists(std::unordered_map<int, std::pair<std::string, int>>& kmers,
                            std::vector<int> kmer,
                            std::unordered_map<int, std::string>& num2str,
                            KMER_DECORATOR kmer_decorator) {
  int hash = get_hash(kmer);
  auto map_entry = kmers.find(hash);
  if(map_entry == kmers.end()) {
    kmers[hash] = std::make_pair(create_kmer(kmer, num2str, kmer_decorator), 0);
  }
}

void update_kmers_with_alphabet(std::unordered_map<int, std::pair<std::string, int>>& kmers,
                                const std::vector<int>& alphabet,
                                std::vector<int>& currentKmer,
                                int k,
                                std::unordered_map<int, std::string>& num2str,
                                KMER_DECORATOR kmer_decorator) {
  if(currentKmer.size() == k) {
    add_kmer_if_not_exists(kmers, currentKmer, num2str, kmer_decorator);
  } else {
    for(char letter: alphabet) {
      currentKmer.push_back(letter);
      update_kmers_with_alphabet(kmers, alphabet, currentKmer, k, num2str, kmer_decorator);
      currentKmer.pop_back();
    }
  }
}

bool is_kmer_allowed(const std::vector<int>& s,
                     const Rcpp::IntegerVector& d,
                     int begin_index,
                     std::unordered_map<int, bool>& isItemAllowed) {
  int current_index = begin_index;
  int i = 0;
  do {
    if(!isItemAllowed[s[current_index]]) {
      return false;
    }
    current_index += d[i] + 1;
    ++i;
  } while (i <= d.length());
  return true;
}

std::unordered_map<std::string, int> __count_kmers(const std::vector<int>& s,
                                                   const Rcpp::IntegerVector& d,
                                                   const std::vector<int>& alphabet,
                                                   std::unordered_map<int, std::string> num2str,
                                                   KMER_DECORATOR kmer_decorator,
                                                   bool pos) {
  std::unordered_map<int, bool> isItemAllowed;
  for(const int& c: alphabet) {
    isItemAllowed[c] = true;
  }
  
  std::unordered_map<int, std::pair<std::string, int>> kmers; // hash -> (kmer, count)
  
  // fill map with found (allowed) kmers
  int window_length = get_window_length(d);
  int last_window_index = s.size() - window_length;
  for(int kmer_begin_index = 0; kmer_begin_index <= last_window_index; ++kmer_begin_index) {
    if(is_kmer_allowed(s, d, kmer_begin_index, isItemAllowed)) {
      update_kmers(kmers, d, s, get_hash(s, d, kmer_begin_index, pos), kmer_begin_index, num2str, kmer_decorator);
    }
  }
  
  if(!pos) {
    // fill map with kmers that can be created with the given alphabet
    std::vector<int> currentKmer;
    currentKmer.clear();
    update_kmers_with_alphabet(kmers, alphabet, currentKmer, (int)(d.size() + 1), num2str, kmer_decorator);
  }
  
  std::unordered_map<std::string, int> res;
  for(const auto& entry: kmers) {
    res[entry.second.first] = entry.second.second;
  }
  return res;
}

template <class B, class S>
void fill_items_coding_maps(B& elems,
                            std::unordered_map<S, int>& val2num,
                            std::unordered_map<int, std::string>& num2str,
                            int& lowest_not_used_num,
                            std::function<std::string(S)> val2str_converter) {
  for(int i = 0; i < elems.size(); ++i) {
    S s_item = static_cast<S>(elems[i]);
    if(val2num.find(s_item) == val2num.end()) {
      val2num[s_item] = lowest_not_used_num;
      num2str[lowest_not_used_num] = val2str_converter(s_item);
      ++lowest_not_used_num;
    }
  }
}

template <class SEQ_TYPE, class ELEM_TYPE>
void fill_encoded_int_vector(SEQ_TYPE str_v,
                             std::vector<int>& res,
                             std::unordered_map<ELEM_TYPE, int>& val2int) {
  for(int i = 0; i < str_v.size(); ++i) {
    res[i] = val2int[static_cast<ELEM_TYPE>(str_v[i])];
  }
}

template <class B, class S>
std::unordered_map<std::string, int> get_kmers(B& s,
               Rcpp::IntegerVector& d,
               B& alphabet,
               std::function<std::string(S)> val2str_converter,
               KMER_DECORATOR kmer_decorator,
               bool pos) {
  // create S -> int (and vice versa) coding maps in order to deal with numbers
  int lowest_not_used_num = 1;
  std::unordered_map<S, int> val2num;
  std::unordered_map<int, std::string> num2str;
  fill_items_coding_maps(s, val2num, num2str, lowest_not_used_num, val2str_converter);
  fill_items_coding_maps(alphabet, val2num, num2str, lowest_not_used_num, val2str_converter);
  
  std::vector<int> int_s(s.size());
  fill_encoded_int_vector<B, S>(s, int_s, val2num);
  
  std::vector<int> int_alphabet(alphabet.size());
  fill_encoded_int_vector<B, S>(alphabet, int_alphabet, val2num);
  
  return __count_kmers(int_s, d, int_alphabet, num2str, kmer_decorator, pos);
}

KMER_DECORATOR get_kmer_decorator(bool pos) {
  return pos ?
    [](std::string& s, int p) { return std::to_string(p) + "_" + s; } :
    [](std::string& s, int p) { return s; };
}

bool is_first_true(Rcpp::LogicalVector& v) {
  return static_cast<bool>(v[0]);
}

// [[Rcpp::export]]
std::unordered_map<std::string, int> countt_kmers_str(Rcpp::StringVector& s,
                      Rcpp::IntegerVector& d, 
                      Rcpp::StringVector& alphabet,
                      Rcpp::LogicalVector& pos) {
  bool positional = is_first_true(pos);
  return get_kmers<Rcpp::StringVector, std::string>(s, d, alphabet,
                                             [](std::string s) { return s; },
                                             get_kmer_decorator(positional),
                                             positional);
}

// [[Rcpp::export]]
std::unordered_map<std::string, int> countt_kmers_num(Rcpp::NumericVector& s,
                      Rcpp::IntegerVector& d,
                      Rcpp::NumericVector& alphabet,
                      Rcpp::LogicalVector& pos) {
  bool positional = is_first_true(pos);
  return get_kmers<Rcpp::NumericVector, double>(s, d, alphabet,
                                         [](double d) { return std::to_string(d); },
                                         get_kmer_decorator(positional),
                                         positional);
}
