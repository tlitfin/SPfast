#include "protein2.h"
#include <stdint.h>
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
	vector<string> Tlist, Qlist;
	string fali, idir, odir, fdb;
	int bpairlist=0, bcheck=0;
	vector<string> folds;
}
static string runtype = "";
bool inbin = false;
using namespace PARAMS2;
int DEBUG = 0;

void rdlist2(string sdir, string flist, vector<string> &slist, const string &suffix){
	string line; vector<string> ss;
	ifstream fs(flist.c_str());
	if(sdir.size()>0 && sdir[sdir.size()-1]!='/') sdir += '/';
	while(getline(fs, line)){
		if(line[0] == '#') continue;
		int n = str2dat(line, ss);
		if(n < 1) continue;
		slist.push_back(sdir + ss[0] + suffix);
	}
	fs.close();
}
inline string getdir2(string sdir0){
	string sdir = sdir0;
	if(sdir.size()>0 && sdir[sdir.size()-1]!='/') sdir += '/';
	return sdir;
}

void rdparams2(int argc, char *argv[]){
	string usage = "Usage: RUN"
	;
	if(argc < 2) die(usage.c_str());
	Tlist.clear(); Qlist.clear();
	fali = idir = odir = "";
	int i0 = 1;
	string suffix;
	while(i0 < argc){
		if(i0==1 && argv[i0][0] != '-'){
			Qlist.push_back(argv[i0 ++]);
			Tlist.push_back(argv[i0 ++]);
			continue;
		}
		string opt1 = argv[i0];
        if(opt1 == "-qlist"){
			suffix = "";
			if(argc>i0+3 && argv[i0+3][0]!='-') suffix = argv[i0+3];
			rdlist2(argv[i0+1], argv[i0+2], Qlist, suffix);
			if(suffix != "") i0 ++;
			i0 += 2;
		} else if(opt1 == "-q"){
			for(; i0<argc-1; i0++){
				if(argv[i0+1][0] == '-') break;
				Qlist.push_back(argv[i0+1]);
			}
		}
		else if(opt1 == "-tdb") {fdb = argv[++i0]; runtype="db";}
		else if(opt1 == "-bin") {inbin=true;}
		else die("unknown option: %s", argv[i0]);
		i0 ++;
	}
}

void run(){
    //FILE *fp = stdout;
    Protein2 *p1;
    for(int j=0; j<Qlist.size(); j++){
        string bn = basename(Qlist[j].c_str());
		string fn = Qlist[j] + ".bin";
		//if (file_existed(fn)) continue;
		p1 = new Protein2();
        p1->rdideal(Qlist[j]);
		p1->wrbin(fn);
        delete p1;
    }
}

void rundb(bool bin){
    //FILE *fp = stdout;
    Protein2 *p1;
	string name;
	std::ofstream ofs(fdb, std::ios::binary);
	int n_entries = Qlist.size();
	ofs.write(reinterpret_cast<const char*>(&n_entries), sizeof(int));
	//int offset = sizeof(int), new_offset=0;
	int64_t offset = (int64_t)sizeof(int), new_offset=(int64_t)0;
    for(int j=0; j<n_entries; j++){
        string bn = basename(Qlist[j].c_str());
		//string fn = Qlist[j] + ".bin";
		//if (file_existed(fn)) continue;
		name = bn + "\n";
		ofs.write(name.c_str(), name.size());
		p1 = new Protein2();
        if (bin){
            p1->rdbin(Qlist[j]);
        } else {
            p1->rdideal(Qlist[j]);
        }
		new_offset = p1->wrbin1(ofs);
		//printf("%s %d\n", bn.c_str(), offset);
		printf("%s %" PRId64 "\n", bn.c_str(), offset);
		offset+=(new_offset+name.size());
        delete p1;
		//fflush(stdout);
		//printf("%s %d %d %d\n", bn.c_str(), offset, new_offset, name.size());
    }
}

int main(int argc, char *argv[]){
	rdparams2(argc, argv);
	if(Qlist.size() < 1) die("no query selected");
	//fflush(stdout);
    if(runtype == "") {
		run();
	} else if (runtype == "db"){
		rundb(inbin);
	}
}
