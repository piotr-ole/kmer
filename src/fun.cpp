// [[Rcpp::plugins("cpp11")]]

#include<Rcpp.h>
#include<string>
#include<unordered_map>
#include<vector>
#include<iostream>
#include<algorithm>

const int HASH_CONST = 101;

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
             int begin_index) {
  long long lres = s[begin_index];
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
                        std::unordered_map<int, std::string>& num2str) {
  int total_kmer_size = get_total_size_of_kmer(s, d, begin_index, num2str);
  std::string res;
  res.reserve(total_kmer_size);
  res += num2str[s[begin_index]];
  int current_index = begin_index;
  for(int d_i = 0; d_i < d.size(); ++d_i) {
    current_index += d[d_i] + 1;
    res += "." + num2str[s[current_index]];
  }
  return res;
}

std::string create_kmer(const std::vector<int>& kmer,
                        std::unordered_map<int, std::string>& num2str) {
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
  return res;
}

void update_kmers(std::unordered_map<int, std::pair<std::string, int>>& kmers,
                  const Rcpp::IntegerVector& d,
                  const std::vector<int>& s,
                  int kmer_hash,
                  int kmer_begin_index,
                  std::unordered_map<int, std::string>& num2str) {
  auto map_elem = kmers.find(kmer_hash);
  if(map_elem == kmers.end()) {
    kmers[kmer_hash] = std::make_pair(create_kmer(s, d, kmer_begin_index, num2str), 0);
  }
  ++kmers[kmer_hash].second;
}

void add_kmer_if_not_exists(std::unordered_map<int, std::pair<std::string, int>>& kmers,
                            std::vector<int> kmer,
                            std::unordered_map<int, std::string>& num2str) {
  int hash = get_hash(kmer);
  auto map_entry = kmers.find(hash);
  if(map_entry == kmers.end()) {
    kmers[hash] = std::make_pair(create_kmer(kmer, num2str), 0);
  }
}

void update_kmers_with_alphabet(std::unordered_map<int, std::pair<std::string, int>>& kmers,
                                const std::vector<int>& alphabet,
                                std::vector<int>& currentKmer,
                                int k,
                                std::unordered_map<int, std::string>& num2str) {
  if(currentKmer.size() == k) {
    add_kmer_if_not_exists(kmers, currentKmer, num2str);
  } else {
    for(char letter: alphabet) {
      currentKmer.push_back(letter);
      update_kmers_with_alphabet(kmers, alphabet, currentKmer, k, num2str);
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
                                                   std::unordered_map<int, std::string> num2str) {
  std::unordered_map<int, bool> isItemAllowed;
  for(const int& c: s) {
    isItemAllowed[c] = true;
  }
  
  std::unordered_map<int, std::pair<std::string, int>> kmers; // hash -> (kmer, count)
  
  // fill map with found (allowed) kmers
  int window_length = get_window_length(d);
  int last_window_index = s.size() - window_length;
  for(int kmer_begin_index = 0; kmer_begin_index <= last_window_index; ++kmer_begin_index) {
    if(is_kmer_allowed(s, d, kmer_begin_index, isItemAllowed)) {
      update_kmers(kmers, d, s, get_hash(s, d, kmer_begin_index), kmer_begin_index, num2str);
    }
  }
  
  // fill map with kmers that can be created with the given alphabet
  std::vector<int> currentKmer;
  currentKmer.clear();
  
  update_kmers_with_alphabet(kmers, alphabet, currentKmer, (int)(d.size() + 1), num2str);
  
  std::unordered_map<std::string, int> res;
  for(const auto& entry: kmers) {
    res[entry.second.first] = entry.second.second;
  }
  return res;
}

void fill_items_coding_maps(Rcpp::StringVector& elems,
                            std::unordered_map<std::string, int>& str2int,
                            std::unordered_map<int, std::string>& num2str,
                            int& lowest_not_used_num) {
  for(int i = 0; i < elems.size(); ++i) {
    std::string str_item = Rcpp::as<std::string>(elems[i]);
    if(str2int.find(str_item) == str2int.end()) {
      str2int[str_item] = lowest_not_used_num;
      num2str[lowest_not_used_num] = str_item;
      ++lowest_not_used_num;
    }
  }
}

void fill_encoded_int_vector(std::vector<std::string> str_v,
                             std::vector<int>& res,
                             std::unordered_map<std::string, int>& str2int) {
  int i = 0;
  for(const std::string& elem: str_v) {
    res[i++] = str2int[elem];
  }
}

// [[Rcpp::export]]
void countt_kmers(Rcpp::StringVector s,
                  Rcpp::IntegerVector d, 
                  Rcpp::StringVector alphabet) {
  // create str -> int (and vice versa) coding maps in order to deal with numbers
  int lowest_not_used_num = 1;
  std::unordered_map<std::string, int> str2int;
  std::unordered_map<int, std::string> num2str;
  fill_items_coding_maps(s, str2int, num2str, lowest_not_used_num);
  fill_items_coding_maps(alphabet, str2int, num2str, lowest_not_used_num);
  
  std::vector<std::string> v_s = Rcpp::as<std::vector<std::string>>(s);
  std::vector<int> int_s(v_s.size());
  fill_encoded_int_vector(v_s, int_s, str2int);

  std::vector<std::string> v_alphabet = Rcpp::as<std::vector<std::string>>(alphabet);
  std::vector<int> int_alphabet(v_alphabet.size());
  fill_encoded_int_vector(v_alphabet, int_alphabet, str2int);
  
  std::unordered_map<std::string, int> res = __count_kmers(int_s, d, int_alphabet, num2str);
  for(const auto& p: res) {
    std::cout << p.first << " " << p.second << std::endl;
  }
}


