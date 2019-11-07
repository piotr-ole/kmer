// [[Rcpp::plugins("cpp11")]]

#include <Rcpp.h>
#include<string>
#include<unordered_map>
#include<vector>

const int HASH_CONST = 101;

const int MOD = 1e9 + 33;

int get_window_length(const Rcpp::IntegerVector& d) {
  int res = 0;
  for(const int& item: d) {
    res += item + 1;
  }
  return res;
}

int get_hash(const std::string& s,
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

std::string create_kmer(const std::string& s,
                        const Rcpp::IntegerVector& d,
                        int begin_index) {
  std::string res;
  int kmer_length = d.size() + 1;
  res.reserve(kmer_length);
  res.append(1, s[begin_index]);
  int current_letter_index = begin_index;
  for(int d_i = 0; d_i < d.size(); ++d_i) {
    current_letter_index += d[d_i] + 1;
    res.append(1, s[current_letter_index]);
  }
  return res;
}

void update_kmers(std::unordered_map<int, std::pair<std::string, int>>& kmers,
                  const Rcpp::IntegerVector& d,
                  const std::string& s,
                  int kmer_hash,
                  int kmer_begin_index) {
  auto map_elem = kmers.find(kmer_hash);
  if(map_elem == kmers.end()) {
    kmers[kmer_hash] = std::make_pair(create_kmer(s, d, kmer_begin_index), 0);
  }
  ++kmers[kmer_hash].second;
}

void update_kmers_with_alphabet(std::unordered_map<int, std::pair<std::string, int>>& kmers,
                                const std::string& alphabet,
                                std::vector<char>& currentKmer,
                                int k) {
  if(currentKmer.size() == k) {
    long long hash = 0;
    for(int i = 0; i < k; ++i) {
      hash = (hash * HASH_CONST) + alphabet[i];
    }
    auto map_entry = kmers.find(hash);
    if(map_entry == kmers.end()) {
      kmers[hash] = std::make_pair(std::string(currentKmer.begin(), currentKmer.end()), 0);
    }
  }
  
  for(char letter: alphabet) {
    currentKmer.push_back(letter);
    update_kmers_with_alphabet(kmers, alphabet, currentKmer, k);
    currentKmer.pop_back();
  }
}

std::unordered_map<std::string, int> count_kmers(Rcpp::CharacterVector s,
                                                 Rcpp::IntegerVector d,
                                                 Rcpp::CharacterVector alphabet) {
  // TODO: filter kmers that have chars not included in the alphabet
  std::string str = Rcpp::as<std::string>(s);
  std::unordered_map<int, std::pair<std::string, int>> kmers; // hash -> (kmer, count)
  int window_length = get_window_length(d);
  int last_window_index = s.size() - window_length;
  for(int kmer_begin_index = 0; kmer_begin_index <= last_window_index; ++kmer_begin_index) {
    update_kmers(kmers, d, str, get_hash(str, d, kmer_begin_index), kmer_begin_index);
  }
  
  std::string str_alphabet = Rcpp::as<std::string>(alphabet);
  std::vector<char> currentKmer;
  update_kmers_with_alphabet(kmers,
                             str_alphabet,
                             currentKmer,
                             (int)(d.size() + 1));
  
  std::unordered_map<std::string, int> res;
  for(const auto& entry: kmers) {
    res[entry.second.first] = entry.second.second;
  }
  return res;
}


// [[Rcpp::export]]
double fun() {
  return 5;
}
