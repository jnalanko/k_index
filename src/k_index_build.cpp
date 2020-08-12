#include <iostream>
#include "r_index.hpp"
#include "build_document_start_marks.hh"
#include <sdsl/rank_support.hpp>
#include <sdsl/rmq_succinct_sada.hpp>
#include <string>
#include <fstream>
#include <streambuf>

using namespace std;
using namespace sdsl;
typedef long long LL;

string read_file_to_string(string filename){

    std::ifstream t(filename);
    std::string str;

    t.seekg(0, std::ios::end);   
    str.reserve(t.tellg());
    t.seekg(0, std::ios::beg);

    str.assign((std::istreambuf_iterator<char>(t)),
                std::istreambuf_iterator<char>());
    return str;
}

LL bit_size(LL x){
    assert(x >= 0);
    LL size = 0;
    while(x != 0){
        x >>= 1;
        size++;
    }
    return size;
}

// LCP of T[iA..iA+k-1] and T[iB..iB+k-1]
LL get_lcp(int_vector<8>& T, LL iA, LL iB, LL k){
    LL L = 0;
    while(L < k && max(iA,iB) + L < T.size() && T[iA + L] == T[iB + L]) L++;
    return L;
}

pair<sd_vector<>, int_vector<> > build_zones_and_LCP(sdsl::int_vector_buffer<>& sa, int_vector<8>& T, LL k){
    LL n = sa.size();
    assert(T.size() == n);

    bit_vector zones(n); // Marks starting positions of k-mer zones
    int_vector<> LCP(1,0,bit_size(k)); // LCP[i] = lcp of k-mers i and i-1 in lex-order. LCP[0] = 0. Initialized to size 1 so buffer doubling works

    LL SA_now = -1;
    LL SA_prev = -1; // Store the previous SA value so we can do a fully sequential scan of the suffix array
    LL kmer_id = 0; // Id of next k-mer to be added
    for(LL i = 0; i < n; i++){
        SA_prev = SA_now;
        SA_now = sa[i];
        // Populate k-mer zones and LCP
        LL lcp = (i == 0) ? 0 : get_lcp(T, SA_now, SA_prev, k);
        if(lcp < k){
            zones[i] = 1;
            if(kmer_id >= LCP.size()) LCP.resize(LCP.size() * 2); // Make room
            LCP[kmer_id] = lcp;
            kmer_id++;
        }
    }
    LCP.resize(kmer_id);

    sd_vector<> zones_sparse(zones);
    return {zones_sparse, LCP};
}

pair<int_vector<>, bit_vector> build_SAS(sdsl::int_vector_buffer<>& sa, sd_vector<>& zones, int_vector<8>& T, LL k){
    LL n = sa.size();
    assert(T.size() == n);
    
    int_vector<> SAS(1,0,bit_size(n)); // Suffix array samples of the last position of every k-mer zone that is not a run-end
    bit_vector marks(1); // SAS_marks[i] = 1 iff k-mer zone i does not end in a run end.
    // ^ These two are initialized to size 1 so that size doubling to make space works

    LL marks_idx = 0;
    LL SAS_idx = 0;
    LL SA_now = -1;
    LL SA_next = -1;
    for(LL i = 0; i < n; i++){

        // Read SA values. We will never back track in the SA so disk streaming works fine.
        // (We read values in order sa[0], sa[1], sa[1], sa[2], sa[2], sa[3]...)
        SA_now = sa[i];
        if(i != n-1) SA_next = sa[i+1];

        // Are we at a run end?
        bool run_end = true;
        if(i != n-1){
            char bwt_i = T[(SA_now - 1 + n) % n];
            char bwt_i_plus_1 = T[(SA_next - 1 + n) % n];
            if(bwt_i == bwt_i_plus_1) run_end = false;
        }

        // Are we at a zone end?
        bool zone_end = (i == n-1 || zones[i+1] == 1);

        // Push to SAS if appropriate
        if(zone_end && !run_end){
            if(SAS_idx >= SAS.size()) SAS.resize(SAS.size() * 2);
            SAS[SAS_idx] = SA_now;
            SAS_idx++;
        }

        if(zone_end){
            // Mark whether we pushed to SAS
            if(marks_idx >= marks.size()) marks.resize(marks.size() * 2);
            marks[marks_idx] = (zone_end && !run_end);
            marks_idx++;
        }
    }

    SAS.resize(SAS_idx);
    marks.resize(marks_idx);
    return {SAS, marks};
}

// T is the text
void SA_callback(sdsl::int_vector_buffer<>& sa, int_vector<8>& T, LL k, string out_prefix){

    cerr << "Building k-mer components" << endl;

    LL n = sa.size();
    assert(T.size() == n);

    sd_vector<> zones; // Marks starting positions of k-mer zones
    int_vector<> LCP; // LCP[i] = lcp of k-mers i and i-1 in lex-order. LCP[0] = 0.
    std::tie(zones,LCP) = build_zones_and_LCP(sa,T,k);

    int_vector<> SAS; // Suffix array samples of the last position of every k-mer zone that is not a run-end
    bit_vector SAS_marks; // SAS_marks[i] = 1 iff k-mer zone i does not end in a run end.
    std::tie(SAS,SAS_marks) = build_SAS(sa,zones,T,k);

    rank_support_v<> SAS_marks_rs(&SAS_marks);
    rmq_succinct_sada<> RMQ(&LCP);

    sd_vector doc_start_marks = build_document_start_marks(T);

    sdsl::store_to_file(zones, out_prefix + ".kmers.zones");
    sdsl::store_to_file(LCP, out_prefix + ".kmers.lcp");
    sdsl::store_to_file(SAS, out_prefix + ".kmers.sas");
    sdsl::store_to_file(SAS_marks, out_prefix + ".kmers.sas.marks");
    sdsl::store_to_file(SAS_marks_rs, out_prefix + ".kmers.sas.marks.rs");
    sdsl::store_to_file(RMQ, out_prefix + ".kmers.rmq"); // todo: rename to lcp.rmq
    sdsl::store_to_file(doc_start_marks, out_prefix + ".docs.starts");

}

int main(int argc, char** argv){

    LL k = 30; // TODO: TO FILE
    string infile = argv[1];
    string data = read_file_to_string(infile);

    auto callback_wrapper = [&k,&infile](sdsl::int_vector_buffer<>& sa, int_vector<8>& T){
        return SA_callback(sa,T,k,infile);
    };

    ri::r_index<> RI(data, true, callback_wrapper);
    ofstream ri_out(infile + ".ri");
    RI.serialize(ri_out);

}