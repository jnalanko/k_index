#pragma once
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/rmq_succinct_sada.hpp>
#include "throwing_streams.hh"
#include "input_reading.hh"
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>

typedef long long LL;
using namespace std;
using namespace sdsl;

// Also clears the file
void check_writable(string filename){
    throwing_ofstream F(filename, std::ofstream::out | std::ofstream::app); // Throws on failure
}

sd_vector<> build_document_start_marks(int_vector<8>& data){
    
    LL n = data.size();
    bit_vector marks(n,0);
    for(LL i = 0; i < n; i++){
        if(i == 0 || (i != n-1 && data[i-1] == '$')) marks[i] = 1;
    }

    sd_vector marks_compressed(marks);
    return marks_compressed;

}
