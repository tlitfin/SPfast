#ifndef _PROTEIN
#define _PROTEIN

#include "sp_misc.h"
#include "sp_type.h"
#include <unordered_set>
#include <stdint.h>

struct hashFunction
{
  size_t operator()(const tuple<int, int, int, int> x) const
  {
    return get<0>(x) ^ get<1>(x) ^ get<2>(x) ^ get<3>(x);
  } 
};

class Protein2{
    string name;
    int nres;
    int nseg;
    vector<int> seq0, resid;
    vector<string> resnum;
    vector<Xvec> x;
    vector<Xvec> x_ideal;
    vector<int> ssec;
    vector<int> segstart;
    vector<int> segmid;
    vector<int> segend;
    vector<int> seglen;
public:
    string info_err;
    Protein2(){};
    Protein2(string fn){rdbin(fn);};
    Protein2(int nr, string seq1, vector<vector<double> > xp);
    Protein2(int nr, string seq1){};
    void rdpdb(string fn);
    void rdideal(string fn);
    void wrpdb(string fn);
    void wrbin(string fn);
    //void wrbin1(std::ofstream &ofs);
    //int wrbin1(std::ofstream &ofs);
    int64_t wrbin1(std::ofstream &ofs);
    void rdbin(string fn);
    void rdbin1(std::ifstream &ifs);
    void wrpdb(FILE*);
    int getnres(){return nres;}
    double distance2(int,int);
    void calSS();
    void calSS2(vector<int> &ssec);
    void add_ideal(int nr, vector<vector<double> > xp);
    void add_dummy_ideal();
    void calineib(double r2);
    void calnneib0(vector<int> &idx);
    void calnneib(vector<int> &idx, vector<int> &aligned, double);
    string getname(){return name;}
    void setname(string fn){name=fn;}
    vector<int> getseq0(){return seq0;}
    vector<Xvec> getx(){return x;}
    vector<int> getssec(){return ssec;}
    vector<int> getsegstart(){return segstart;}
    vector<int> getsegmid(){return segmid;}
    vector<int> getsegend(){return segend;}
    ~Protein2(){};
    friend class Salign;
};
void quatfit(int n, vector<double> &w, vector<Xvec> &x1, vector<Xvec> &x2, double u[][3]);
void translate(int na, vector<Xvec> &x0, vector<Xvec> &xn, double *xc);
void translate2(int na, vector<double> &wfit, vector<Xvec> &x0, vector<Xvec> &xn, double *xc);
void rotmol(int na, double u[][3], vector<Xvec> &x0, vector<Xvec> &xn);
void rotmol_ideal(int nsa, double u[][3], vector<Xvec> &x0, vector<Xvec> &xn, vector<int> &start, vector<int> &mid, vector<int> &end, vector<int> &ssec);
inline double distance2(double *xap, double *xbp);
//
//
class Salign{
    int ma=0, mb=0, msa=0, msb=0;
    Protein2 *pa, *pb;
    bool breverse;
    int nA, nB, nsA, nsB, Lmin, seeds=0, valid_seeds=0;
    vector<Xvec> *xap, *xbp;        // quote for abbrev.
    vector<Xvec> *xap_ideal, *xbp_ideal;        // quote for abbrev.
//
    double D0, D0_2;                // Rcut^2 from len_ali
    vector<Xvec> xn1, xn2;      // for changed coords
    vector<double> wfit;
    double coarse=0.;
// for DP
    double gap0, gap1;
    double ss_prefilter;
    double Rscore, u[4][3];     //  current score; rotation matrix
    double scores_all[20] = {0.};       // series of scores saved as ieLA... below
    TwoDVector<double> Rmat, Smat, Rmat_fragment, Smat_fragment;
    TwoDVector<short> Idir, Idir_fragment;
    vector<int> ialign[3];      // the aligned idx, only 2 are used 
// saved & max
    double Rscore_sv, u_sv[4][3];
    vector<int> ialign_sv[3];
    double Rscore_max, u_max[4][3];
    vector<int> ialign_max[3];
    vector<int> ialign_fragment[3];
    vector<int> ialign_fragment_sv[3];
    vector<int> ialign_fragment_max[3];
    vector<int> list2[2];
    unordered_set<tuple<int, int, int, int>, hashFunction> success3;
    int imax_bounded, jmax_bounded;
    double dmax_bounded;
    bool filtered = false;
    double xc1[3], xc2[3];
public:
    Salign():ma(0),mb(0){}
    ~Salign(){delete_x();}
    int init(Protein2 *a, Protein2 *b, int brev);
    void init(Protein2 *a, Protein2*b){init(a,b,0);};
    double run_pairwise(Protein2 *a, Protein2*b);
    double get_coarse(){return coarse;};
    void initalign(string smat[], string sali[]);
    void delete_x();
    void alloc_x();
    void alloc_lengths();
    void alloc_large();
    double optimize_score(int,int);
    double optimize_align();
    void init_block_align();
    double optimize_block_align();
    void Run_all();
    void Run1();
    void Run2();
    void RunSS(double mismatch, double matcha=1., double matchb=1.);
    void RunSS_fragment(double mismatch, double matcha=1., double matchb=1.);
    double PrepareSeed(int i, int j, int offset_i, int offset_, double rms0);
    void RunSeg();
    void RunSeg3();
    void Run3();
    void Run4();
    double run_Align1();
    double run_Align1_fragment(int anchor[][3]);
    double Align_DP();
    double Align_DP_bounded(int i_start, int i_end, int j_start, int j_end);
    void compute_forward(int i, int j, double &dmax, int &imax, int &jmax);
    void score_segment(int i, int j);
    double Align_DP_fragment(bool score_only);
    double Align_DP_fragment_bounded(int anchor[][3]);
    void traceback_DP_bounded();
    void block_align(vector<int> ialign_cache[3]);
    void refine_align(vector<int> ialign_cache[3]);
    int calrms(double &rms, double cut2);
    int redo_fit(int flag, double c);
    int redo_fit_ideal(int flag, double c);
    double redo_fit_ideal_seed(double minScore);
    int fit1(vector<int> *iali2);
    void calRmat();
    void calRmat_fragment();
    void calRmat_fragment_bounded(int anchor[][3]);
    void calRmat_bounded(int i0, int i1, int j0, int j1);
    void calRmat_ss();
    double calRscore();
    double calRscore_frag();
    double getscore(){return Rscore;}
    void calscores_all();
    void calLali2(int*, double);
    bool save_align();
    void score_final();
    void print_max(FILE*);
    void print_max(string&);
    void restore_max();
    void prtali(FILE*);
    void prtdis_pair(string &sinfo);
    void prtali(string&);
    void checkAli1(string);
    void checkAlign(string idir, string fn, FILE *fpo, vector<Protein2*>&, Protein2 *p0);
    void pro_swap(int);
    void scoreOnly();
    void calSPscores(double*);
    double calSPscore(){double SP[4]; calSPscores(SP); return SP[0];};
    double calSPscoreQ(){double SP[4]; calSPscores(SP); return SP[1];};
    double getscoresall(){return scores_all[1];};
    double getRmax(){return Rscore_max;};
    bool isFiltered(){return filtered;}
    double (*get_u())[4][3]{ double (*ptr)[4][3]; ptr = &u; return ptr;}
};
enum {iSP, iTM, iLG, iGDT};
enum {ieSP, ieTMa, ieLG, ieGDT, ieLA, ieRMS, ieLE, ieTMb, ieTMc, ieSEQ};
//
//                 TMc, TMa, TMb if iTM
namespace Salign_PARAMS{
    extern int iprint, inorm, score_type, fragsize, bscoreOnly;
    extern double Alpha, Beta, D00, Denv, cutoff, gap_min, scaling_factor;
    extern double segcut, segalpha, ssprefcut, pref_gap0, pref_gap1, final_gap0, convergence_criterion, coarsecut;
    extern int riters;
    extern bool bfullalign, bsingledom;
    extern string outpdb;
}
using namespace Salign_PARAMS;
//
//
static const int nres_std = 20;
static const char rnam3_std[][4] = 
        {"ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU",
        "MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"};
static const char rnam1_std0[] = "XACDEFGHIKLMNPQRSTVWY";
static const char *rnam1_std = rnam1_std0 + 1;
inline int aaDefine(string rn0, int DEBUG=0){
    string rn = rn0;
    if(rn=="HSD" || rn=="HSE") rn = "HIS";
    for(int i=0; i<nres_std; i++){
        if(rn==rnam3_std[i]) return i;
    }
    if(DEBUG > 0){
        fprintf(stderr, "Unrecognized residue name: %s\n", rn.c_str());
    }
    return -1;
}
inline int aaDefine1(char rn, int DEBUG=0){
    for(int i=0; i<nres_std; i++){
        if(rn==rnam1_std[i]) return i;
    }
    if(DEBUG > 0){
        fprintf(stderr, "Unrecognized residue name: %c\n", rn);
    }
    return -1;
}

#endif
