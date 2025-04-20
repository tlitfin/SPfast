#include "protein2.h"
#include <inttypes.h>

namespace Salign_PARAMS2{
// alpha is the normalized factor(as seen in paper); D00 is "D0"; Denv the distance cutoff used to determine the "environment residues"
    double Alpha=0.3, D00=4., Denv=12., cutoff=-1, gap_min=0., scaling_factor=3.75;
// score_type decids the type of calculated alignment scores (SP, TM, GDT-scores are supported)
// fragsize is the length of fragment for the original alignment trials during finding the best alignment
    int score_type=iSP, fragsize=20; 
// iprint controls the details to print, the bigger the more printed
    int iprint=-9999;
// bscoreOnly: 0, structure alignment; others, score only (using predefined alignment)
// 1,scored according to resi No.; 2,scored according to residues sequentially
    int bscoreOnly=0;
    bool bfullalign=0;
}
namespace PARAMS2{
    using namespace Salign_PARAMS2;
//  vector<string> Tlist, Qlist;
    unordered_set<string> Qset;
    string fali, idir, odir, fdb, fdbout;
    int bpairlist=0, bcheck=0;
    vector<string> folds;
}
static string runtype = "";
using namespace PARAMS2;
int DEBUG = 0;

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
inline string getdir2(string sdir0){
    string sdir = sdir0;
    if(sdir.size()>0 && sdir[sdir.size()-1]!='/') sdir += '/';
    return sdir;
}

void rdparams2(int argc, char *argv[]){
    string usage = "Usage: RUN fdb odir [-q query / -qlist query_list] [-tdb in.db]";
    if(argc < 2) die(usage.c_str());
    fali = idir = odir = "";
    int i0 = 1;
    string suffix;
    while(i0 < argc){
        if(i0==1 && argv[i0][0] != '-'){
            fdb = argv[i0++];
            odir = argv[i0++];
            continue;
        }
        string opt1 = argv[i0];
        if(opt1 == "-qlist"){
            suffix = "";
            if(argc>i0+2 && argv[i0+2][0]!='-') suffix = argv[i0+2];
            rdlist2(argv[i0+1], Qset, suffix);
            if(suffix != "") i0 ++;
            i0 ++;
        } else if(opt1 == "-q"){
            for(; i0<argc-1; i0++){
                if(argv[i0+1][0] == '-') break;
                Qset.insert(argv[i0+1]);
            }
        }
		else if(opt1 == "-tdb") {fdbout = argv[++i0]; runtype="db";}
        else die("unknown option: %s", argv[i0]);
        i0 ++;
    }
}

void run(){
    //FILE *fp = stdout;
    std::ifstream ifs(fdb, std::ios::binary);
    std::ifstream ifs_index(fdb + ".index", std::ios::binary);
    if(ifs.fail()) die("Invalid DB");
    if(ifs_index.fail()) die("Invalid DB index");

    Protein2 Tpro;
    string out_fn, name, line;
    vector<string> ss;
    //char* endptr;
    while(getline(ifs_index, line)){
        int n = str2dat(line, ss);
        //printf("%s\n", ss[0].c_str());
        if(Qset.find(ss[0]) != Qset.end()){
            //printf("%s %d\n", ss[0].c_str(), stoi(ss[1]));
            Tpro = Protein2();
            //ifs.seekg(stoi(ss[1]));
            //ifs.seekg(strtoimax(ss[1], &endptr, 10));
            ifs.seekg(stoll(ss[1]));
            std::getline(ifs, name);
            Tpro.setname(name);
            Tpro.rdbin1(ifs);
            out_fn = odir + "/"+ name + ".bin";
            Tpro.wrbin(out_fn);
            Qset.erase(ss[0]);
            if (Qset.size()==0) break;
        }
    }
}

void rundb(){
    //indb
	std::ifstream ifs(fdb, std::ios::binary);
    std::ifstream ifs_index(fdb + ".index", std::ios::binary);
    if(ifs.fail()) die("Invalid DB");
    if(ifs_index.fail()) die("Invalid DB index");

	//outdb
    std::ofstream ofs(fdbout, std::ios::binary);
    int n_entries = Qset.size();
    ofs.write(reinterpret_cast<const char*>(&n_entries), sizeof(int));
    
	Protein2 *p1;
    string out_fn, name, line, bn;
    vector<string> ss;

    int64_t offset = (int64_t)sizeof(int), new_offset=(int64_t)0;
    while(getline(ifs_index, line)){
		int n = str2dat(line, ss);
        if(Qset.find(ss[0]) != Qset.end()){
			p1 = new Protein2();

			ifs.seekg(stoll(ss[1]));
			std::getline(ifs, bn);
			name = bn + "\n";
			
			ofs.write(name.c_str(), name.size());
			p1->rdbin1(ifs);

			new_offset = p1->wrbin1(ofs);
			printf("%s %" PRId64 "\n", bn.c_str(), offset);
			offset+=(new_offset+name.size());
			delete p1;
    	}
	}
}

int main(int argc, char *argv[]){
    rdparams2(argc, argv);
    if(runtype == "") {
        run();
    } else if (runtype == "db"){
		rundb();
	}
}
