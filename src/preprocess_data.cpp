#include "input_reading.hh"

#include <string>
#include <iostream>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;
using namespace std;

int64_t replaced_count = 0;

void upper_case_and_replace_non_ACGT(string& S){
    for(int64_t i = 0; i < S.size(); i++){
        char c = toupper(S[i]);
        if(c != 'A' && c != 'C' && c != 'G' && c != 'T'){
            int64_t r = rand() % 4;
            if(r == 0) c = 'A';
            if(r == 1) c = 'C';
            if(r == 2) c = 'G';
            if(r == 3) c = 'T';
            replaced_count++;
        }
        S[i] = c;
    }
}

int main(int argc, char** argv){

    if(argc != 3){
        cerr << "Usage: " << argv[0] <<  "genomes_path outfile" << endl;
        return 1;
    }

    string genomes_path = argv[1];
    string outfile = argv[2];

    throwing_ofstream concat_out(outfile);
    throwing_ofstream doc_names_out(outfile + ".docs.names");

    for (const auto & entry : fs::directory_iterator(genomes_path)){
        string path = entry.path().string();
        doc_names_out << path << "\n";
        cerr << "Processing " << path << endl;

        Sequence_Reader sr(path, FASTA_MODE);
        string genome;
        while(!sr.done()){
            genome += sr.get_next_query_stream().get_all();
        }
        upper_case_and_replace_non_ACGT(genome);
        concat_out << genome << "$";
    }
//    out << "#";

    cerr << "Replaced " << replaced_count << " non-ACGT chars" << endl;

}
