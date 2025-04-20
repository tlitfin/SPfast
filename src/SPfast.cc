#include "protein2.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <thread>
#include <unordered_map>
#include <mutex>
#include <condition_variable>
#include <unistd.h>

namespace Salign_PARAMS{
// alpha is the normalized factor(as seen in paper); D00 is "D0"; Denv the distance cutoff used to determine the "environment residues"
//  double Alpha=0.3, D00=4., Denv=12., cutoff=-1, gap_min=0., scaling_factor=3.75; //SP0
    double Alpha=0.4, D00=5.5, Denv=12., cutoff=-1, gap_min=0., scaling_factor=5.625; //SPf
//TOM params
    double segcut=5.0, segalpha=0.5, ssprefcut=-1., pref_gap0=1.0, pref_gap1=0.5, final_gap0=0.2, convergence_criterion=0.05, coarsecut=-1.;
    //double ssprefcut=0.35; //SCOPe benchmark
    int riters=1;
// score_type decids the type of calculated alignment scores (SP, TM, GDT-scores are supported)
// fragsize is the length of fragment for the original alignment trials during finding the best alignment
    int score_type=iSP, fragsize=0; 
// iprint controls the details to print, the bigger the more printed
    int iprint=-9999;
// bscoreOnly: 0, structure alignment; others, score only (using predefined alignment)
// 1,scored according to resi No.; 2,scored according to residues sequentially
    int bscoreOnly=0;
    bool bfullalign=0;
    bool bsingledom=0;
    int batchstart=0, batchend=-1;
}
namespace PARAMS{
    using namespace Salign_PARAMS;
    vector<string> Tlist, Qlist;
    unordered_set<string> Tset;
    string fali, idir, odir, fdb;
    int bpairlist=0, bcheck=0;
    vector<string> folds;
}
static string runtype = "";
using namespace PARAMS;
int DEBUG = 0;
void rdlist(string sdir, string flist, vector<string> &slist, const string &suffix){
    string line; vector<string> ss;
    ifstream fs(flist.c_str());
    if(sdir.size()>0 && sdir[sdir.size()-1]!='/') sdir += '/';
    while(getline(fs, line)){ //gcc/6.3.0
        if(line[0] == '#') continue;
        int n = str2dat(line, ss);
        if(n < 1) continue;
        slist.push_back(sdir + ss[0] + suffix);
    }
    fs.close();
}
void rdlist2(string flist, unordered_set<string> &slist, const string &suffix){
    string line; vector<string> ss;
    ifstream fs(flist.c_str());
    while(getline(fs, line)){
        if(line[0] == '#') continue;
        int n = str2dat(line, ss);
        if(n < 1) continue;
        slist.insert(ss[0]+suffix);
    }
    fs.close();
}
inline string getdir(string sdir0){
    string sdir = sdir0;
    if(sdir.size()>0 && sdir[sdir.size()-1]!='/') sdir += '/';
    return sdir;
}
void rdparams(int argc, char *argv[]){
    string usage = "Usage: RUN [-pair pdb1 pdb2|-pairlist sdir list] [-iprint 1] [-SP|-TM|-GDT] [-fast] [-cutoff cut]"
    ;
    if(argc < 3) die(usage.c_str());
    Tlist.clear(); Qlist.clear();
    fali = idir = odir = "";
    int i0 = 1;
    string suffix;

    while(i0 < argc){
        if(i0==1 && argv[i0][0] != '-'){
            Qlist.push_back(argv[i0 ++]);
            Tlist.push_back(argv[i0 ++]);
            iprint = 100;
            continue;
        }
        string opt1 = argv[i0];
        if(opt1 == "-tlist"){   // [-tlist tdir list [suffix]
            suffix = "";
            if(argc>i0+3 && argv[i0+3][0]!='-') suffix = argv[i0+3];
            rdlist(argv[i0+1], argv[i0+2], Tlist, suffix);
            if(suffix != "") i0 ++;
            i0 += 2;
        }else if(opt1 == "-qlist"){
            suffix = "";
            if(argc>i0+3 && argv[i0+3][0]!='-') suffix = argv[i0+3];
            rdlist(argv[i0+1], argv[i0+2], Qlist, suffix);
            if(suffix != "") i0 ++;
            i0 += 2;
        }else if(opt1 == "-subdb"){
            suffix = "";
            if(argc>i0+2 && argv[i0+2][0]!='-') suffix = argv[i0+2];
            rdlist2(argv[i0+1], Tset, suffix);
            if(suffix != "") i0 ++;
            i0 ++;
        }else if(opt1 == "-pairlist"){
            bpairlist = 1;
            suffix = "";
            if(argc>i0+3 && argv[i0+3][0]!='-') suffix = argv[i0+3];
            rdlist(argv[i0+1], argv[i0+2], Tlist, suffix);
            if(suffix != "") i0 ++;
            i0 += 2;
        }else if(opt1 == "-pair"){
            if(iprint == -9999) iprint = 100;
            Qlist.push_back(idir + argv[++ i0]);
            Tlist.push_back(idir + argv[++ i0]);
        }else if(opt1 == "-t"){
            for(; i0<argc-1; i0++){
                if(argv[i0+1][0] == '-') break;
                Tlist.push_back(argv[i0+1]);
            }
        } else if(opt1 == "-q"){
            for(; i0<argc-1; i0++){
                if(argv[i0+1][0] == '-') break;
                Qlist.push_back(argv[i0+1]);
            }
        }
        else if(opt1 == "-idir") idir = getdir(argv[++i0]);
        else if(opt1 == "-odir") odir = getdir(argv[++i0]);
        else if (opt1 == "-iprint") iprint = atoi(argv[++i0]);
        else if(opt1 == "-scoreOnly") bscoreOnly = 1;
        else if(opt1 == "-scoreOnly2") bscoreOnly = 2;
        else if (opt1 == "-a") Alpha = strtod(argv[++i0], NULL);
        else if (opt1 == "-d0") D00 = strtod(argv[++i0], NULL);
        else if (opt1 == "-segcut") segcut = strtod(argv[++i0], NULL);
        else if (opt1 == "-segalpha") segalpha = strtod(argv[++i0], NULL);
        else if (opt1 == "-ssprefcut") ssprefcut = strtod(argv[++i0], NULL);
        else if (opt1 == "-coarsecut") coarsecut = strtod(argv[++i0], NULL);
        else if (opt1 == "-prefgap0") pref_gap0 = strtod(argv[++i0], NULL);
        else if (opt1 == "-finalgap0") final_gap0 = strtod(argv[++i0], NULL);
        else if (opt1 == "-converge") convergence_criterion = strtod(argv[++i0], NULL);
        else if (opt1 == "-prefgap1") pref_gap1 = strtod(argv[++i0], NULL);
        else if (opt1 == "-singledom") bsingledom = 1;
        else if (opt1 == "-riters") riters = atoi(argv[++i0]);
        else if (opt1 == "-gapmin") gap_min = strtod(argv[++i0], NULL);
        else if(opt1 == "-frag") fragsize = atoi(argv[++i0]);
        else if(opt1 == "-fast") {convergence_criterion=0.9;coarsecut=5.5;segcut=4.0;}
        else if(opt1 =="-SPscore") {Alpha=0.3; D00=4.; scaling_factor=3.75;}
        else if (opt1 == "-denv") Denv = strtod(argv[++i0], NULL);
//      else if (opt1 == "-outpdb") outpdb = argv[++i0];
        else if(opt1 == "-check") {fali = argv[++i0]; bcheck=1;}
        else if(opt1 == "-fali") {fali = argv[++i0]; bcheck=2;}
        else if(opt1 == "-fadir") {fali = argv[++i0]; bcheck=3;}
        else if(opt1 == "-cutoff") cutoff = strtod(argv[++i0], NULL);
        else if(opt1 == "-plist") {fali = argv[++i0]; runtype="plist";}
        else if(opt1 == "-tdb") {fdb = argv[++i0]; runtype="db";}
        else if(opt1 == "-SP") score_type = iSP;
        else if(opt1 == "-GDT") score_type = iGDT;
        else if(opt1 == "-TM") score_type = iTM;
        else if(opt1 == "-LG") score_type = iLG;
        else if(opt1 == "-fullalign") bfullalign = 1;
        else if(opt1 == "-batchstart") batchstart = atoi(argv[++i0]);
        else if(opt1 == "-batchend") batchend = atoi(argv[++i0]);
        else die("unknown option: %s", argv[i0]);
        i0 ++;
    }
    if(iprint == -9999) iprint = 0;
    if(bpairlist) Qlist = Tlist;
    //fprintf(stderr, "protein number: %d %d\n", Qlist.size(), Tlist.size());
}

void align_query_parallel(vector<Protein2> &Tpro){
    int n0 = 0, ntpl=Tpro.size();
    // scan the list of "query" proteins
    #if defined MP
        #pragma omp parallel for schedule(dynamic, 1)
        //#pragma omp parallel for schedule(static, 1)
    #endif
    for(int j=0; j<Qlist.size(); j++){
        FILE *fp;
        Protein2 *p1;
        string bn = basename(Qlist[j].c_str());
        p1 = new Protein2(Qlist[j]);
        if(odir != ""){
            string fn = odir + bn + ".sp";
            fp = openfile(fn, "w");
        } else {
            fp = stdout;
        }

        string* sdim = new string[ntpl]; //2.3M was too big for stack (before new)- segfault - may not be neccessary at all - must optimize!!
        for(int i=0; i<ntpl; i++) sdim[i] = "";
        for(int i=n0; i<ntpl; i++){
            Salign *sa1 = new Salign();
            int err = sa1 -> init(p1, &Tpro[i], 1);
            if(! err) {
            sa1 -> Run_all();
            if (sa1 -> isFiltered()){delete sa1; continue;}
#if defined MP
            sa1 -> print_max(sdim[i]); //if MP and not odir
#else
            sa1 -> print_max(fp);
#endif
            }
            delete sa1;
        }
#if defined MP
        for(int i=n0; i<ntpl; i++){
            fprintf(fp, "%s", sdim[i].c_str());
        }
#endif

        if(odir != "") fclose(fp); //b
        if(! bpairlist) delete p1;
        delete[] sdim;
    }
}

void align_tpl_parallel(vector<Protein2> &Tpro){
    int n0 = 0, ntpl=Tpro.size();

    int chunk_size = max(1, ntpl/2000);

    // scan the list of "query" proteins
    for(int j=0; j<Qlist.size(); j++){
        FILE *fp;
        Protein2 *p1;
        string bn = basename(Qlist[j].c_str());
        p1 = new Protein2(Qlist[j]);
        if(odir != ""){
            string fn = odir + bn + ".sp";
            fp = openfile(fn, "w");
        } else {
            fp = stdout;
        }
        //printf("%s\n", p1->getname().c_str());

        string* sdim = new string[ntpl]; //2.3M was too big for stack (before new)- segfault - may not be neccessary at all - must optimize!!
        for(int i=0; i<ntpl; i++) sdim[i] = "";
        #if defined MP
            #pragma omp parallel for schedule(dynamic, chunk_size)
        #endif
        for(int i=n0; i<ntpl; i++){
            Salign *sa1 = new Salign();
            int err = sa1 -> init(p1, &Tpro[i], 1);
            if(! err) {
            sa1 -> Run_all();
            if (sa1 -> isFiltered()){delete sa1; continue;}
#if defined MP
            sa1 -> print_max(sdim[i]); //if MP and not odir
#else
            sa1 -> print_max(fp);
#endif
            }
            delete sa1;
        }
#if defined MP
        for(int i=n0; i<ntpl; i++){
            fprintf(fp, "%s", sdim[i].c_str());
        }
#endif

        if(odir != "") fclose(fp); //b
        if(! bpairlist) delete p1;
        delete[] sdim;
    }
}

void align_list_mp(){
// read the list of "template" proteins
    vector<Protein2> Tpro(Tlist.size());

    int n0 = 0, ntpl=Tlist.size(); 
    #if defined MP
        #pragma omp parallel for schedule(dynamic, 1)
    #endif
    for(int i=0; i<ntpl; i++){
        Tpro[i] = Protein2(Tlist[i]);
    }

    unsigned int nthreads = std::thread::hardware_concurrency();
    if (Qlist.size()>=(5*nthreads)){
        align_query_parallel(Tpro);
    } else {
        align_tpl_parallel(Tpro);
    }   
}

void align_db(string fdb, int batchstart=0, int batchend=-1){
    std::ifstream ifs(fdb, std::ios::binary);
    if(ifs.fail()) die("Invalid DB");
    std::ifstream ifs_index(fdb + ".index", std::ios::binary);
    if(ifs_index.fail()) die("Invalid DB index");

    int npro_db;
    ifs.read(reinterpret_cast<char*>(&npro_db), sizeof(npro_db));
    vector<string> ss;
    string line, name;
    vector<Protein2> Tpro;
    int ntpl, count;

    //If searching subdb
    if (Tset.size()>0){
        ntpl = Tset.size();
        count = 0;
        Tpro.resize(ntpl);
 
        while(getline(ifs_index, line)){
            int n = str2dat(line, ss);
            if(Tset.find(ss[0]) != Tset.end()){
                ifs.seekg(stoll(ss[1]));
                std::getline(ifs, name);
                Tpro[count].setname(name);
                Tpro[count].rdbin1(ifs);
                count++;
            }
        }
    } else {
        //if batching
        if ((batchstart>0) && (batchend>0)){
            //find starting offset
            count = -1;
            while(getline(ifs_index, line)){
                count+=1;
                if (count>=batchstart){
                    int n = str2dat(line, ss);
                    ifs.seekg(stoll(ss[1]));
                    break;
                }
            }
        }

        // read the list of "template" proteins
        if (batchend>0)
            if (batchend<=npro_db){
                ntpl=batchend-batchstart;
            } else {
                ntpl=npro_db-batchstart;
            }
        else {
            ntpl=npro_db;
        }

        Tpro.resize(ntpl);      
        for(int i=0; i<ntpl; i++){
            Tpro[i] = Protein2();
            //sequential read - no seek required
            std::getline(ifs, name);
            Tpro[i].setname(name);
            Tpro[i].rdbin1(ifs);
        }
    }

    unsigned int nthreads = std::thread::hardware_concurrency();
    if (Qlist.size()>=(5*nthreads)){
        align_query_parallel(Tpro);
    } else {
        align_tpl_parallel(Tpro);
    }   
}

void read_pairs(unordered_set<string> &proteins, unordered_map<string, Protein2*> &Plist, std::mutex &mtx, std::condition_variable &cv){
    for (auto p1=proteins.begin(); p1!=proteins.end(); ++p1){
        string fn1 = idir + *p1 + ".ideal.bin";
        std::lock_guard<std::mutex> lock(mtx);
        if(Plist.count(fn1) < 1){
            Plist[fn1] = new Protein2(fn1);
            cv.notify_all();
        }
    }
}

Protein2 *getData(string fn, unordered_map<string, Protein2*> &Plist, std::mutex &mtx, std::condition_variable &cv){
    std::unique_lock<std::mutex> lock(mtx); //it will block here until an iteration of readfile is complete independent of if it needs to block
    cv.wait(lock, [&fn, &Plist]() {
        return Plist.count(fn) > 0; //condition is checked first - wont wait/block if already true
    });
    return Plist[fn];
}

// only run for file containing pairs
void align_plist(string idir, string fali){
    unordered_map<string, Protein2*> Plist;
    std::mutex mtx;
    std::condition_variable cv;
    
    unordered_set<string> proteins;
    vector<tuple<string, string>> pairs;
    FILE *fp = fopen(fali.c_str(), "r");
    char str[501], p1[501], p2[501];
    while(fgets(str, 500, fp) != NULL){
        sscanf(str, "%s%s", p1, p2);
        proteins.insert(p1);
        proteins.insert(p2);
        pairs.push_back(tuple<string, string>{p1, p2});
    }
    fclose(fp);

    std::thread readerThread(read_pairs, std::ref(proteins), std::ref(Plist), std::ref(mtx), std::ref(cv));

    #if defined MP
        #pragma omp parallel for
    #endif
    for (int i=0; i<pairs.size(); i++){
        Salign *sa1 = new Salign();

        string p1 = std::get<0>(pairs[i]);
        string p2 = std::get<1>(pairs[i]);
        string fn1 = idir + p1 + ".ideal.bin", fn2 = idir + p2 + ".ideal.bin";
        
        sa1 -> init(getData(fn1, Plist, mtx, cv), getData(fn2, Plist, mtx, cv), 1);
        sa1 -> Run_all();
        if (sa1 -> isFiltered()) {delete sa1; continue;}
        sa1 -> print_max(stdout);
        delete sa1;
    }
    readerThread.join();
}
int main(int argc, char *argv[]){
    ios_base::sync_with_stdio(false);
    rdparams(argc, argv);
    if(Tlist.size() < 1 && runtype == "") die("no tpl selected");
    if(Qlist.size() < 1 && runtype == "") die("no query selected");
    if(runtype == "") {
        align_list_mp();
    } else if(runtype == "plist"){
        align_plist(idir, fali);
    } else if(runtype == "db"){
        if ((batchend>0) && (batchend<batchstart)) die("Invalid batch");
        align_db(fdb, batchstart, batchend);
    } else die("unknown running mode: %s", runtype.c_str());
}
