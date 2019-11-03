#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cstdlib>
#include <string>
#include <map>
#include <list>
#include <utility>
#include <exception>

using namespace Rcpp;


// [[Rcpp::export]]
std::list<std::string> map_up_to_k(std::string sequence, unsigned int k)
{
  std::list<std::string> l = std::list<std::string>();
  for (unsigned int len = 1; len <= k; len++) 
  {
    for (unsigned int i = 0; i <= sequence.size() - len; i++)
    {
      l.push_back(sequence.substr(i, len));
    }
  }
  return l;
}

// [[Rcpp::export]]
std::list<std::string> map(std::string sequence, unsigned int k)
{
  std::list<std::string> l = std::list<std::string>();
  for (unsigned int i = 0; i <= sequence.size() - k; i++)
  {
    l.push_back(sequence.substr(i, k));
  }
  return l;
}

// [[Rcpp::export]]
std::map<std::string, int> reduce(std::list<std::string> l) {
  std::map<std::string, int> mapa = std::map<std::string, int>();
  std::list<std::string>::iterator it = l.begin();
  for (it = l.begin() ; it != l.end(); it++)
  {
    try
    {
      mapa.at((*it))++;
    }
    catch (std::out_of_range& e)
    {
      mapa.insert(std::pair<std::string, int>((*it) , 1));
    }
  }
  return mapa;
}


// [[Rcpp::export]]
std::map<std::string, int> map_reduce(std::string sequence, int k)
  {
  std::map<std::string, int> mapa = std::map<std::string, int>();
  for (unsigned int i = 0; i <= sequence.size() - k; i++)
  {
    std::string s = sequence.substr(i, k);
    try
    {
      mapa.at(s)++;
    }
    catch (std::out_of_range& e)
    {
      mapa.insert(std::pair<std::string, int>(s , 1));
    }
  }
  return mapa;
  }
