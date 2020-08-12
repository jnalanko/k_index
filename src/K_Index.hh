#pragma once

#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/rmq_succinct_sada.hpp>
#include "throwing_streams.hh"
#include "input_reading.hh"
#include "r_index.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>

using namespace std;
typedef long long LL;

class K_Index{

private: // No copying

  K_Index(const K_Index&) = delete;
  void operator=(const K_Index&) = delete;

public:

    ri::r_index<> RI;
    sd_vector<> zones; // Marks starting positions of k-mer zones
    rank_support_sd<> zones_rs;
    select_support_sd<> zones_ss;
    int_vector<> LCP; // LCP between every distinct k-mers adjacent in the lex-order
    rmq_succinct_sada<> RMQ; // RMQ for the k-mer LCP
    int_vector<> SAS; // Suffix array samples of the last position of every k-mer zone that is not a run-end
    bit_vector SAS_marks; // SAS_marks[i] = 1 iff k-mer zone i does not end in a run end.
    rank_support_v<> SAS_marks_rs;
    sd_vector<> doc_starts; // Marks the first text position of every document in the concatenation
    rank_support_sd<> doc_starts_rs;
    select_support_sd<> doc_starts_ss;
    LL k;

    K_Index(LL k) : k(k) {}

    void load(string filename_prefix){
        cerr << "Loading r-index" << endl;

	    std::ifstream in(filename_prefix + ".ri");

        RI.load(in);

        cerr << "Loading kmer zones..." << endl;
        
        sdsl::load_from_file(zones, filename_prefix + ".kmers.zones");
        sdsl::util::init_support(zones_rs, &zones);
        sdsl::util::init_support(zones_ss, &zones);

        cerr << "Loading kmer LCP..." << endl;
        
        sdsl::load_from_file(LCP, filename_prefix + ".kmers.lcp");

        cerr << "Loading kmer LCP RMQ..." << endl;
        
        sdsl::load_from_file(RMQ, filename_prefix + ".kmers.rmq");

        cerr << "Loading k-mer SA samples..." << endl;

        sdsl::load_from_file(SAS, filename_prefix + ".kmers.sas");

        cerr << "Loading k-mer SA samples marks..." << endl;

        sdsl::load_from_file(SAS_marks, filename_prefix + ".kmers.sas.marks");
        sdsl::load_from_file(SAS_marks_rs, filename_prefix + ".kmers.sas.marks.rs");
        SAS_marks_rs.set_vector(&SAS_marks);

        cerr << "Loading document start marks..." << endl;

        sdsl::load_from_file(doc_starts, filename_prefix + ".docs.starts");
        sdsl::util::init_support(doc_starts_rs, &doc_starts);
        sdsl::util::init_support(doc_starts_ss, &doc_starts);

        cerr << "Index loading finished" << endl;
    }

    // id is between 0 and (number of kmers - 1)
    pair<LL,LL> kmer_id_to_kmer_zone(LL id){
        LL n_kmers = LCP.size();
        LL start = zones_ss(id + 1);
        LL end = id == n_kmers - 1 ? RI.bwt_size() - 1 : zones_ss(id+2) - 1;
        return {start, end};
    }

    // Returns the 0-based id the k-mer zone containing lex
    LL lex_pos_to_kmer_id(LL lex){
        return zones_rs.rank(lex+1) - 1; // K-mer id of current k-mer
    }

    // Check if the k-mer can be left extended by c such that c is the only possible
    // left-extension and the resulting k-mer has the same frequency as the current one.
    bool is_redundant(LL kmer_id, char c){
        // A k-mer is redundant if all characters in its k-mer zone are the
        // same and left-extending leads to a k-mer zone of the same size.

        LL n_kmers = LCP.size();

        LL zone_start, zone_end;
        std::tie(zone_start, zone_end) = kmer_id_to_kmer_zone(kmer_id);

        if(RI[zone_start] != c) 
            return false; // Wrong extension

        LL c_count = RI.bwt_rank(zone_end + 1, c) - RI.bwt_rank(zone_start, c);

        if(c_count != zone_end - zone_start + 1) 
            return false; // Two or more left extensions

        LL p = RI.LF(zone_start);
        LL pred_kmer_id = lex_pos_to_kmer_id(p);
        LL zone_start_pred, zone_end_pred;
        std::tie(zone_start_pred, zone_end_pred) = kmer_id_to_kmer_zone(pred_kmer_id);
        return zone_end - zone_start == zone_end_pred - zone_start_pred;
    }

    // pos needs to be a run end or a k-mer zone end
    LL get_SA_sample(LL pos){
        if(pos == RI.bwt_size() - 1 || RI[pos] != RI[pos+1]){
            // At run end
            return RI.get_next_SA_sample(pos);
        } else{
            // At k-mer zone end
            assert(zones[pos+1] == 1);
            LL kmer_id = lex_pos_to_kmer_id(pos);
            LL sample_id = SAS_marks_rs.rank(kmer_id);
            assert(sample_id < SAS.size());
            return SAS[sample_id];
        }
    }

    // Takes MS_k[i] and the lex-position of one match of length at least MK_k[i]
    // Returns MS_k[i-1] and a lex-position of one match of length at least MS_k[i-1]
    // By updating the given parameters match_length and lex
    void MS_k_update(LL c, LL& match_length, LL& lex){
        if(RI[lex] == c){ // Left extension success
            lex = RI.LF(lex);
            match_length = min(k,match_length+1);
        } else{ // Left extension failure
            LL p1 = RI.bwt_pred(lex,c); // Returns -1 if does not exist
            assert(p1 < lex);
            LL p2 = RI.bwt_succ(lex,c); // Returns bwt_size if does not exist
            assert(p2 > lex);
            if(p1 == -1 && p2 == RI.bwt_size()){
                // c is not found in the index
                lex = 0;
                match_length = 0;
                return;
            }

            LL q_now = lex_pos_to_kmer_id(lex);

            LL lcp1 = -1;
            LL lcp2 = -1;
            if(p1 != -1){ // Predecessor exists
                LL q1 = lex_pos_to_kmer_id(p1);
                if(q1 == q_now) lcp1 = k; // In the same k-mer
                else lcp1 = LCP[RMQ(q1+1, q_now)];
            }
            if(p2 != RI.bwt_size()){ // Successor exists
                LL q2 = lex_pos_to_kmer_id(p2);
                if(q2 == q_now) lcp2 = k; // In the same k-mer
                else lcp2 = LCP[RMQ(q_now+1, q2)];
            }

            match_length = min(match_length, max(lcp1, lcp2));

            if(lcp1 > lcp2) lex = p1;
            else lex = p2;

            lex = RI.LF(lex);
            match_length = min(k,match_length+1);
        }
    }

    vector<LL> locate_kmer_occs(LL kmer_id){
        vector<LL> occs;
        if(kmer_id == -1) return occs;

        LL start_of_range, end_of_range;
        std::tie(start_of_range, end_of_range) = kmer_id_to_kmer_zone(kmer_id);
        LL n_occs = end_of_range - start_of_range + 1;

        LL text_pos = get_SA_sample(end_of_range);
        for(LL j = 0; j < n_occs; j++){
            occs.push_back(text_pos);
            text_pos = RI.Phi(text_pos);
        }

        return occs;

    }

    // Returns a vector v where v[i] is the id of the k-mer starting from read[i], or -1
    // if the k-mer is not found in the index or i >= read.size()-k+1
    vector<LL> search_kmers(string& read){
        vector<LL> ids;
        LL lex = 0; // Position of an occurrence
        LL match_length = 0;
        for(LL i = read.size()-1; i >= 0; i--){
            
            MS_k_update(read[i], match_length, lex); // Updates the values of match_length and lex

            LL kmer_id = lex_pos_to_kmer_id(lex); // k-mer id of current k-mer

            if(match_length == k) ids.push_back(kmer_id); // k-mer found
            else ids.push_back(-1); // k-mer not found
        }

        std::reverse(ids.begin(), ids.end());
        return ids;
    }

    // Returns pair (document id, position in document)
    pair<LL, LL> locate_in_document(LL text_pos){
        LL doc_id = doc_starts_rs(text_pos+1) - 1; // -1: to 0-based indexing
        LL doc_start = doc_starts_ss(doc_id+1);
        return {doc_id, text_pos - doc_start};
    }

};

