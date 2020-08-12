#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/rmq_succinct_sada.hpp>
#include "throwing_streams.hh"
#include "K_Index.hh"
#include "input_reading.hh"
#include "r_index.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>

using namespace std;
typedef long long LL;

string figure_out_file_format(string filename){
    for(LL i = filename.size()-1; i >= 0; i--){
        if(filename[i] == '.'){
            string end = filename.substr(i);
            
            if(end == ".fasta") return "fasta";
            if(end == ".fna") return "fasta";
            if(end == ".ffn") return "fasta";
            if(end == ".faa") return "fasta";
            if(end == ".frn") return "fasta";

            if(end == ".fastq") return "fastq";
            if(end == ".fq") return "fastq";

            if(end == ".gz") return "gzip";

            throw(runtime_error("Unknown file format: " + filename));
        }
    }
    throw(runtime_error("Unknown file format: " + filename));
    return "unknown";
}

// Stores the intersection into buf1 and returns the number of elements in the
// intersection (does not resize buf1). Buffer elements must be sorted.
// Assumes all elements in a buffer are distinct
LL intersect_buffers(vector<LL>& buf1, LL buf1_len, vector<LL>& buf2, LL buf2_len){

    LL i = 0, j = 0, k = 0;
    while(i < buf1_len && j < buf2_len){
        if(buf1[i] < buf2[j]) i++;
        else if(buf1[i] > buf2[j]) j++;
        else{
            buf1[k] = buf1[i];
            i++; j++; k++;
        }
    }
    return k;

}


// Returns in sorted order
vector<LL> hits_to_colorset(K_Index& KI, vector<LL>& hits){

    vector<LL> C;
    for(LL text_pos : hits){
        C.push_back(KI.locate_in_document(text_pos).first);
    }
    std::sort(C.begin(), C.end());
    C.erase(std::unique(C.begin(), C.end()), C.end()); // Delete duplicates
    return C;

}

int main(int argc, char** argv){

    string index_prefix = argv[1];
    string query_file = argv[2];

    LL mode = figure_out_file_format(query_file) == "fasta" ? FASTA_MODE : FASTQ_MODE;
    Sequence_Reader sr(query_file, mode);

    K_Index KI(30); // todo: load k from file
    KI.load(index_prefix);

    LL skipped = 0;
    LL not_skipped = 0;
    LL read_id = 0;
    while(!sr.done()){

        string read = sr.get_next_query_stream().get_all();
        vector<LL> kmer_ids = KI.search_kmers(read);
        vector<LL> intersection;
        LL intersection_size = -1;

        // Identify documents and update intersection
        for(LL i = 0; i < read.size(); i++){
            LL kmer_id = kmer_ids[i];
            if(kmer_id == -1){
                skipped++;
                continue;
            }
            if(i != 0 && KI.is_redundant(kmer_id, read[i-1])){
                skipped++;
                continue; // kmer i has the same color set as kmer i-1
            }
            not_skipped++;
            vector<LL> hits = KI.locate_kmer_occs(kmer_id);
            vector<LL> C = hits_to_colorset(KI, hits);
            if(intersection_size == -1){
                intersection = C;
                intersection_size = C.size();
            }
            else{
                intersection_size = intersect_buffers(intersection, intersection_size, C, C.size());
            }
        }
        
        cout << read_id << " ";
        if(intersection_size != -1){
            for(LL i = 0; i < intersection_size; i++){
                cout << intersection[i] << " ";
            }
        }
        cout << "\n";
        read_id++;
    }

    cerr << "Fraction skipped " << (double)skipped / (skipped + not_skipped) << endl;
    
}






