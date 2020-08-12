#pragma once

#include <cmath>

typedef long long LL;

class Progress_printer{

    public:

    LL n_jobs;
    LL processed;
    LL total_prints;
    LL next_print;

    Progress_printer(LL n_jobs, LL total_prints) : n_jobs(n_jobs), processed(0), total_prints(total_prints), next_print(0) {}

    void job_done(){
        if(next_print == processed){
            LL progress_percent = round(100 * ((double)processed / n_jobs));
            cerr << "Progress: " << progress_percent << + "%" << endl;
            next_print += n_jobs / total_prints;
        }
        processed++;
    }

};

