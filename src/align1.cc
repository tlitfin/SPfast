#include "protein2.h"
#include <unordered_set>
#include "quatfit_theo.h"


const double GDT_r2[] = {1, 2*2, 4*4, 8*8};
inline double getd0_TM(double L0){
    if(L0 > 15) return max(.5, 1.24 * pow(L0-15., 1./3) - 1.8);
    return 0.5;
}
const double Eff0 = 1.3;    // to include those pairs between 8~9A
double Salign::Align_DP(){
    gap0 = max(gap_min, gap0);
    gap1 = max(gap_min, gap1);
    //double match_open_penalty = 0.1;
    
    double xd[3], dat[3];
    for(int j=0; j<nB+1; j++) {Smat[0][j] = 0; Idir[0][j] = -1;}
    for(int i=0; i<nA+1; i++) {Smat[i][0] = 0; Idir[i][0] = -1;}
    Idir[0][0] = 0;
    double dmax = -100.;
    int imax=-1, jmax=-1;
    for(int i=1; i<nA+1; i++)
    for(int j=1; j<nB+1; j++){
        double dt = Rmat[i-1][j-1];
        dat[0] = Smat[i-1][j-1] + dt;
        dat[1] = Smat[i][j-1];
        dat[2] = Smat[i-1][j];
        if(Idir[i][j-1] != 1) dat[1] -= gap0;       // simple DP
        else dat[1] -= gap1;
        if(Idir[i-1][j] != 2) dat[2] -= gap0;
        else dat[2] -= gap1;

        //if(Idir[i][j-1] != 0) dat[0] -= match_open_penalty;   // NEW to discourage isolated matches distorting alignment and inflating Le
        int it = 0;
        if(dat[it] < dat[1]) it = 1;
        if(dat[it] < dat[2]) it = 2;
        Smat[i][j] = dat[it]; Idir[i][j] = it;
        if(dat[it] <= 0) { Smat[i][j] = 0; Idir[i][j] = -1; }
        if(dmax < Smat[i][j]) {dmax=Smat[i][j]; imax=i; jmax=j;}
    }

// trace
    int imax0 = imax, jmax0 =jmax;
    if(imax<0 || jmax<0){
        printf("%d %d %f\n", imax, jmax, dmax);
        die("not used smat: %s %s", pa->name.c_str(), pb->name.c_str());
    }
    for(int m=0; m<2; m++) ialign[m].clear();
    int i0, j0, k0;
    while(imax>=0 || jmax>=0){
        if(imax==0 && jmax==0) break;
        k0 = Idir[imax][jmax];
        i0 = imax - 1; j0 = jmax - 1;
        if(k0 == 0) {imax --; jmax --;}
        else if(k0 == 1) {jmax --; i0 = -1;}
        else if(k0 == 2) {imax --; j0 = -1;}
        else if(k0 == -1) break;
        else die("not known k0: %d\n", k0);
        ialign[0].push_back(i0);
        ialign[1].push_back(j0);
    }
    for(int m=0; m<2; m++) reverse(ialign[m].begin(), ialign[m].end());
    return dmax;
}

void Salign::init_block_align(){
    //TODO: optimize reset - only need to reset bounds of blocks - resetting is >5% by valgrind (now >8%)
    Idir.reset();
    Smat.reset();
    Rmat.reset();

    for(int j=0; j<nB+1; j++) {Smat[0][j] = 0; Idir[0][j] = -1;}
    for(int i=0; i<nA+1; i++) {Smat[i][0] = 0; Idir[i][0] = -1;}
    Idir[0][0] = 0;
    dmax_bounded = -100.;
    imax_bounded=-1, jmax_bounded=-1;
}

double Salign::Align_DP_bounded(int i_start, int i_end, int j_start, int j_end){
    gap0 = max(gap_min, gap0);
    gap1 = max(gap_min, gap1);
    //double match_open_penalty = 0.1;
    double xd[3], dat[3];

    for(int i=(i_start+1); i<=(i_end+1); i++)
    for(int j=(j_start+1); j<=(j_end+1); j++){
        double dt = Rmat[i-1][j-1];
        dat[0] = Smat[i-1][j-1] + dt;
        dat[1] = Smat[i][j-1];
        dat[2] = Smat[i-1][j];
        if(Idir[i][j-1] != 1) dat[1] -= gap0;       // simple DP
        else dat[1] -= gap1;
        if(Idir[i-1][j] != 2) dat[2] -= gap0;
        else dat[2] -= gap1;

        //if(Idir[i][j-1] != 0) dat[0] -= match_open_penalty;       // NEW to discourage isolated matches distorting alignment and inflating Le
        int it = 0;
        if(dat[it] < dat[1]) it = 1;
        if(dat[it] < dat[2]) it = 2;
        Smat[i][j] = dat[it]; Idir[i][j] = it;
        if(dat[it] <= 0) { Smat[i][j] = 0; Idir[i][j] = -1; } //went <=0 for N-term
        if(dmax_bounded < Smat[i][j]) {dmax_bounded=Smat[i][j]; imax_bounded=i; jmax_bounded=j;}
    }
    return dmax_bounded;
}

void Salign::traceback_DP_bounded(){
    // trace
    int imax0 = imax_bounded, jmax0 = jmax_bounded;
    if(imax_bounded<0 || jmax_bounded<0){
        printf("%d %d %f\n", imax_bounded, jmax_bounded, dmax_bounded);
        die("not used smat: %s %s", pa->name.c_str(), pb->name.c_str());
    }
    for(int m=0; m<2; m++) ialign[m].clear();
    int i0, j0, k0;
    while(imax_bounded>=0 || jmax_bounded>=0){
        if(imax_bounded==0 && jmax_bounded==0) break;
        k0 = Idir[imax_bounded][jmax_bounded];
        i0 = imax_bounded - 1; j0 = jmax_bounded - 1;
        if(k0 == 0) {imax_bounded --; jmax_bounded --;}
        else if(k0 == 1) {jmax_bounded --; i0 = -1;}
        else if(k0 == 2) {imax_bounded --; j0 = -1;}
        else if(k0 == -1) break;
        else die("not known k0: %d\n", k0);
        ialign[0].push_back(i0);
        ialign[1].push_back(j0);
    }
    for(int m=0; m<2; m++) reverse(ialign[m].begin(), ialign[m].end());
}

double Salign::Align_DP_fragment(bool score_only){
    auto& startA = pa->segstart;
    auto& midA = pa->segmid;
    auto& endA = pa->segend;
    auto& startB = pb->segstart;
    auto& midB = pb->segmid;
    auto& endB = pb->segend;
    double xd[3], dat[3];
    for(int j=0; j<nsB+1; j++) {Smat_fragment[0][j] = 0; Idir_fragment[0][j] = -1;}
    for(int i=0; i<nsA+1; i++) {Smat_fragment[i][0] = 0; Idir_fragment[i][0] = -1;}
    Idir_fragment[0][0] = 0;
    double dmax = -100.;
    int imax=-1, jmax=-1;
    for(int i=1; i<nsA+1; i++)
    for(int j=1; j<nsB+1; j++){
        double dt = Rmat_fragment[i-1][j-1];
        dat[0] = Smat_fragment[i-1][j-1] + dt; // - 0.3;
        dat[1] = Smat_fragment[i][j-1];
        dat[2] = Smat_fragment[i-1][j];
        if(Idir_fragment[i][j-1] != 1) dat[1] -= gap0;      // simple DP
        else dat[1] -= gap1;
        if(Idir_fragment[i-1][j] != 2) dat[2] -= gap0;
        else dat[2] -= gap1;
//
        int it = 0;
        if(dat[it] < dat[1]) it = 1;
        if(dat[it] < dat[2]) it = 2;
        Smat_fragment[i][j] = dat[it]; Idir_fragment[i][j] = it;

        if(dat[it] <= 0) { Smat_fragment[i][j] = 0; Idir_fragment[i][j] = -1; }
        if(dmax < Smat_fragment[i][j]) {dmax=Smat_fragment[i][j]; imax=i; jmax=j;}
    }

// trace
    if (score_only) return dmax;
    int imax0 = imax, jmax0 =jmax;
    if(imax<0 || jmax<0){
        printf("%d %d %f\n", imax, jmax, dmax);
        die("not used smat: %s %s", pa->name.c_str(), pb->name.c_str());
    }
    for(int m=0; m<2; m++) ialign[m].clear();
    for(int m=0; m<2; m++) ialign_fragment[m].clear();
    int i0, j0, k0;
    while(imax>=0 || jmax>=0){
        if(imax==0 && jmax==0) break;
        k0 = Idir_fragment[imax][jmax];
        i0 = imax - 1; j0 = jmax - 1;
        if(k0 == 0) {imax --; jmax --;}
        else if(k0 == 1) {jmax --; i0 = -1;}
        else if(k0 == 2) {imax --; j0 = -1;}
        else if(k0 == -1) break;
        else die("not known k0: %d\n", k0);
        if ((i0==-1) || (j0==-1)) continue;
        ialign[0].push_back(endA[i0]);
        ialign[0].push_back(midA[i0]);
        ialign[0].push_back(startA[i0]);
        ialign[1].push_back(endB[j0]);
        ialign[1].push_back(midB[j0]);
        ialign[1].push_back(startB[j0]);

        ialign_fragment[0].push_back(i0);
        ialign_fragment[1].push_back(j0);
    }
    for(int m=0; m<2; m++) reverse(ialign[m].begin(), ialign[m].end());
    for(int m=0; m<2; m++) reverse(ialign_fragment[m].begin(), ialign_fragment[m].end());
    
    return dmax;
}

inline void Salign::compute_forward(int i, int j, double &dmax, int &imax, int &jmax){
    double xd[3], dat[3];
    double dt = Rmat_fragment[i-1][j-1];
    dat[0] = Smat_fragment[i-1][j-1] + dt; // - 0.3;
    dat[1] = Smat_fragment[i][j-1];
    dat[2] = Smat_fragment[i-1][j];
    if(Idir_fragment[i][j-1] != 1) dat[1] -= gap0;      // simple DP
    else dat[1] -= gap1;
    if(Idir_fragment[i-1][j] != 2) dat[2] -= gap0;
    else dat[2] -= gap1;
//
    int it = 0;
    if(dat[it] < dat[1]) it = 1;
    if(dat[it] < dat[2]) it = 2;
    Smat_fragment[i][j] = dat[it]; Idir_fragment[i][j] = it;
    if(dat[it] <= 0) { Smat_fragment[i][j] = 0; Idir_fragment[i][j] = -1; }
    if(dmax < Smat_fragment[i][j]) {dmax=Smat_fragment[i][j]; imax=i; jmax=j;}
}

double Salign::Align_DP_fragment_bounded(int anchor[2][3]){
    //NOTE: we are assuming that the segments must all be included in the alignment - may otherwise cause a problem
    auto& startA = pa->segstart;
    auto& midA = pa->segmid;
    auto& endA = pa->segend;
    auto& startB = pb->segstart;
    auto& midB = pb->segmid;
    auto& endB = pb->segend;

    //initialization
    double xd[3], dat[3];
    Smat_fragment.reset();
    Idir_fragment.fill(-1);
    for(int j=0; j<nsB+1; j++) {Smat_fragment[0][j] = 0; Idir_fragment[0][j] = -1;}
    for(int i=0; i<nsA+1; i++) {Smat_fragment[i][0] = 0; Idir_fragment[i][0] = -1;}
    Idir_fragment[0][0] = 0;
    double dmax = -100.;
    int imax=-1, jmax=-1;

    //pre-anchor
    for(int i=1; i<(anchor[0][0]+1); i++)
    for(int j=1; j<(anchor[1][0]+1); j++){  
        compute_forward(i, j, dmax, imax, jmax);
    }

    //anchor
    int i=anchor[0][0]+1, j=anchor[1][0]+1;
    Smat_fragment[i][j] = Smat_fragment[i-1][j-1] + Rmat_fragment[i-1][j-1]; 
    Idir_fragment[i][j] = 0;

    i=anchor[0][1]+1; j=anchor[1][1]+1;
    if ((i==-1) && (j==-1)){
        //pass
    } else if (i==0){ //was -1
        Smat_fragment[anchor[0][0]+1][j] = Smat_fragment[anchor[0][0]+1][j-1] - gap0;
        Idir_fragment[anchor[0][0]+1][j] = 1;
    } else if (j==0){ //was -1
        Smat_fragment[i][anchor[1][0]+1] = Smat_fragment[i-1][anchor[1][0]+1] - gap0;
        Idir_fragment[i][anchor[1][0]+1] = 2;
    } else {
        Smat_fragment[i][j] = Smat_fragment[i-1][j-1] + Rmat_fragment[i-1][j-1]; 
        Idir_fragment[i][j] = 0;
    }

    i=anchor[0][2]+1, j=anchor[1][2]+1;
    Smat_fragment[i][j] = Smat_fragment[i-1][j-1] + Rmat_fragment[i-1][j-1]; 
    Idir_fragment[i][j] = 0;

    //best path must include anchor
    dmax=Smat_fragment[i][j]; imax=i; jmax=j;

    //post-anchor
    for(int i=anchor[0][2]+1; i<nsA+1; i++)
    for(int j=anchor[1][2]+1; j<nsB+1; j++){    
        compute_forward(i, j, dmax, imax, jmax);
    }

    //traceback
    int imax0 = imax, jmax0 =jmax;
    if(imax<0 || jmax<0){
        printf("%d %d %f\n", imax, jmax, dmax);
        die("not used smat: %s %s", pa->name.c_str(), pb->name.c_str());
    }
    for(int m=0; m<2; m++) ialign[m].clear();
    for(int m=0; m<2; m++) ialign_fragment[m].clear();
    int i0, j0, k0;
    while(imax>=0 || jmax>=0){
        if(imax==0 && jmax==0) break;
        k0 = Idir_fragment[imax][jmax];
        i0 = imax - 1; j0 = jmax - 1;
        if(k0 == 0) {imax --; jmax --;}
        else if(k0 == 1) {jmax --; i0 = -1;}
        else if(k0 == 2) {imax --; j0 = -1;}
        else if(k0 == -1) break;
        else die("not known k0: %d\n", k0);
        if ((i0==-1) || (j0==-1)) continue;
        ialign[0].push_back(endA[i0]);
        ialign[0].push_back(midA[i0]);
        ialign[0].push_back(startA[i0]);
        ialign[1].push_back(endB[j0]);
        ialign[1].push_back(midB[j0]);
        ialign[1].push_back(startB[j0]);

        ialign_fragment[0].push_back(i0);
        ialign_fragment[1].push_back(j0);
    }
    for(int m=0; m<2; m++) reverse(ialign[m].begin(), ialign[m].end());
    for(int m=0; m<2; m++) reverse(ialign_fragment[m].begin(), ialign_fragment[m].end());
    return dmax;
}

int Salign::redo_fit(int flag, double rcut2=64){
// flag: 0: all; 1: w=(1+r2/d2); 2:
    //double xc1[3], xc2[3];
    if(flag == 1 && rcut2 < 0.1) rcut2 = D0_2;
    bool bSP = (score_type == iSP);
    bool bGDT = (score_type == iGDT);
    int nfit = 0;
    for(int m=0; m<2; m++) list2[m].clear();
    for(int i=0; i<ialign[0].size(); i++){
        int i1=ialign[0][i], i2=ialign[1][i];
        if(i1<0 || i2<0) continue;
        if(flag == 0) {
            wfit[nfit] = 1;
        }else {
            double r2 = xn1[i1].distance2( (*xbp)[i2] );
            if(flag == 1) {
                if((bSP ||bGDT) && r2>4*D0_2*Eff0) continue;
                double dt = r2 / rcut2;
                wfit[nfit] = 1 / (1 + dt);
//              wfit[nfit] = 1 / pow(1 + dt,2); // test on 03/21/12, <0.1% diff
            } else if(flag == 2){
                if(r2 > rcut2) continue;
                wfit[nfit] = 1.;
            }
        }
// to protect xn1 from being covered by small nfit
        list2[0].push_back(i1);
        list2[1].push_back(i2);
        nfit ++;
    }
    if(nfit < 3){
        if(DEBUG > 0) fprintf(stderr, "too small alignment size: %d %d\n", ialign[0].size(), nfit);
        return 1;
        rotmol(nA, u, (*xap), *xap);
        pa->wrpdb("1.pdb"); pb->wrpdb("2.pdb");
        exit(0);
    }
//
    assert(nfit == list2[0].size());
    for(int i=0; i<nfit; i++){
        int i1 = list2[0].at(i), i2 = list2[1][i];
        for(int m=0; m<3; m++){
            xn1[i][m] = (*xap)[i1][m];
            xn2[i][m] = (*xbp)[i2][m];
        }
    }
    double ds = 0.;
    for(int i=0; i<nfit; i++) ds += wfit[i];
    if(ds > 1.0e-8){
        for(int i=0; i<nfit; i++) wfit[i] /= ds;
    } else {
        for(int i=0; i<nfit; i++) wfit[i] = 1./nfit;
    }
    translate2(nfit, wfit, xn1, xn1, xc1);
    translate2(nfit, wfit, xn2, xn2, xc2);
    //quatfit(nfit, wfit, xn2, xn1, u);

    ///*
    double A[9];
    double rmsd;
    double minScore=-1;
    double E0 = InnerProduct(A, xn1, xn2, nfit, wfit);
    double rot[3][3];

    double bla = 1.;
    FastCalcRMSDAndRotation(rot, A, &rmsd, E0, bla, minScore);
    for (int m=0; m<3; m++){
        for (int n=0; n<3; n++){
            u[n][m] = rot[m][n];
        }
    }
    //*/

    for(int m=0; m<3; m++) u[3][m] = xc2[m] - dot_product(u[m], xc1);
    rotmol(nA, u, (*xap), xn1);
    return 0;
}

//redo_fit_ideal - maintain separate xap(for ideal coords)
int Salign::redo_fit_ideal(int flag, double rcut2=64){
// flag: 0: all; 1: w=(1+r2/d2); 2:
    //double xc1[3], xc2[3];
    if(flag == 1 && rcut2 < 0.1) rcut2 = D0_2;
    bool bSP = (score_type == iSP);
    bool bGDT = (score_type == iGDT);
    int nfit = 0;
    int i1, i2;
    for(int m=0; m<2; m++) list2[m].clear();
    for(int i=0; i<ialign[0].size(); i++){
        i1=ialign[0][i], i2=ialign[1][i];
        if(i1<0 || i2<0) continue;
        if(flag == 0) {
            wfit[nfit] = 1;
        }else {
            double r2 = xn1[i1].distance2( (*xbp_ideal)[i2] );
            if(flag == 1) {
                if((bSP ||bGDT) && r2>4*D0_2*Eff0) continue;
                double dt = r2 / rcut2;
                wfit[nfit] = 1 / (1 + dt);
//              wfit[nfit] = 1 / pow(1 + dt,2); // test on 03/21/12, <0.1% diff
            } else if(flag == 2){
                if(r2 > rcut2) continue;
                wfit[nfit] = 1.;
            }
        }
// to protect xn1 from being covered by small nfit
        list2[0].push_back(i1);
        list2[1].push_back(i2);
        nfit ++;
    }
    if(nfit < 3){
        if(DEBUG > 0) fprintf(stderr, "too small alignment size: %d %d\n", ialign[0].size(), nfit);
        return 1;
        //rotmol_ideal(nsA, u, (*xap_ideal), *xap_ideal, pa->segstart, pa->segmid, pa->segend);
        //pa->wrpdb("1.pdb"); pb->wrpdb("2.pdb");
        //exit(0);
    }
//
    assert(nfit == list2[0].size());
    for(int i=0; i<nfit; i++){
        i1 = list2[0].at(i), i2 = list2[1][i];
        for(int m=0; m<3; m++){
            xn1[i][m] = (*xap_ideal)[i1][m];
            xn2[i][m] = (*xbp_ideal)[i2][m];
        }
    }
    double ds = 0.;
    for(int i=0; i<nfit; i++) ds += wfit[i];
    if(ds > 1.0e-8){
        for(int i=0; i<nfit; i++) wfit[i] /= ds;
    } else {
        for(int i=0; i<nfit; i++) wfit[i] = 1./nfit;
    }
    translate2(nfit, wfit, xn1, xn1, xc1);
    translate2(nfit, wfit, xn2, xn2, xc2);
    //quatfit(nfit, wfit, xn2, xn1, u);

    //
    double A[9];
    double rmsd;
    double minScore=-1;
    double E0 = InnerProduct(A, xn1, xn2, nfit, wfit);
    double rot[3][3];

    double bla = 1.;
    FastCalcRMSDAndRotation(rot, A, &rmsd, E0, bla, minScore);
    
    //u = rot
    for (int m=0; m<3; m++){
        for (int n=0; n<3; n++){
            u[n][m] = rot[m][n];
        }
    }
    
    for(int m=0; m<3; m++) u[3][m] = xc2[m] - dot_product(u[m], xc1);
    rotmol_ideal(nsA, u, (*xap_ideal), xn1, pa->segstart, pa->segmid, pa->segend, pa->ssec);
    return 0;
}

double Salign::redo_fit_ideal_seed(double minScore){
// flag: 0: all; 1: w=(1+r2/d2); 2:
    //double xc1[3], xc2[3];
    int nfit = 0;
    for(int m=0; m<2; m++) list2[m].clear();
    for(int i=0; i<ialign[0].size(); i++){
        int i1=ialign[0][i], i2=ialign[1][i];
        wfit[nfit] = 1;
        list2[0].push_back(i1);
        list2[1].push_back(i2);
        nfit ++;
    }
    if(nfit < 3){
        if(DEBUG > 0) fprintf(stderr, "too small alignment size: %d %d\n", ialign[0].size(), nfit);
        return 1000.;
        //rotmol_ideal(nsA, u, (*xap_ideal), *xap_ideal, pa->segstart, pa->segmid, pa->segend);
        //pa->wrpdb("1.pdb"); pb->wrpdb("2.pdb");
        //exit(0);
    }
//
    assert(nfit == list2[0].size());
    for(int i=0; i<nfit; i++){
        int i1 = list2[0].at(i), i2 = list2[1][i];
        for(int m=0; m<3; m++){
            xn1[i][m] = (*xap_ideal)[i1][m];
            xn2[i][m] = (*xbp_ideal)[i2][m];
        }
    }
    double ds = 0.;
    for(int i=0; i<nfit; i++) ds += wfit[i];
    if(ds > 1.0e-8){
        for(int i=0; i<nfit; i++) wfit[i] /= ds;
    } else {
        for(int i=0; i<nfit; i++) wfit[i] = 1./nfit;
    }
    translate2(nfit, wfit, xn1, xn1, xc1);
    translate2(nfit, wfit, xn2, xn2, xc2);
    //quatfit(nfit, wfit, xn2, xn1, u);

    //
    double A[9];
    double rmsd;
    double E0 = InnerProduct(A, xn1, xn2, nfit, wfit);
    double rot[3][3];

    double bla = 1.;
    FastCalcRMSDAndRotation(rot, A, &rmsd, E0, bla, minScore);
    
    if (rmsd>minScore) return rmsd;

    //u = rot
    for (int m=0; m<3; m++){
        for (int n=0; n<3; n++){
            u[n][m] = rot[m][n];
        }
    }
    for(int m=0; m<3; m++) u[3][m] = xc2[m] - dot_product(u[m], xc1);
    rotmol_ideal(nsA, u, (*xap_ideal), xn1, pa->segstart, pa->segmid, pa->segend, pa->ssec);
    return rmsd;
}

int Salign::calrms(double &rms, double cut2=-1.){
    int n1 = 0; rms = 0.;
    for(int i=0; i<ialign[0].size(); i++){
        int i1=ialign[0][i], i2=ialign[1][i];
        if(i1<0 || i2<0) continue;
        double r2 = xn1[i1].distance2( (*xbp)[i2] );
        if(cut2 > 0 && r2 > cut2) continue;
        rms += r2;
        n1 ++;
    }
    if(n1 > 0) rms = sqrt(rms / n1);
    return n1;
}

// fill Rmat for DP by using "xn1, xbp, D0_2"
void Salign::calRmat(){
    for(int i=0; i<nA; i++)
    for(int j=0; j<nB; j++){
        double r2 = xn1[i].distance2((*xbp)[j]);
        if(score_type == iSP){
            Rmat[i][j] = 0.;
            if(r2 < D0_2*4*Eff0) Rmat[i][j] = 1. / (1. + r2 / D0_2) - 1/(1+4*Eff0);     // slightly extended to include those 8~9A
            //else Rmat[i][j] = -0.3; ////TOM: NEW!!!!!!! to avoid garbage matched residues with long length - should also stop penalizing transition between insert -> delete
        } else if(score_type == iGDT) {
            Rmat[i][j] = 0.;
            if(r2 > GDT_r2[3]){
                if(r2 < GDT_r2[3]*Eff0) Rmat[i][j] = 1.e-3;
                continue;
            }
            for(int m=3; m>=0; m--){
                if(r2 > GDT_r2[m]) break;
                Rmat[i][j] = 1. - m/4.;
            }
        } else {
            Rmat[i][j] = 1. / (1. + r2 / D0_2);
        }
    }
}

void Salign::calRmat_bounded(int i0, int i1, int j0, int j1){
    for(int i=i0; i<=i1; i++){
    for(int j=j0; j<=j1; j++){
        double r2 = xn1[i].distance2((*xbp)[j]);
        if(score_type == iSP){
            Rmat[i][j] = 0.; //otherwise must use this
            //Rmat[i][j] = -1.; ////TOM: NEW!!!!!!! to avoid garbage matched residues with long length - should also stop penalizing transition between insert -> delete (should be relative to gap size)
            if(r2 < D0_2*4*Eff0) Rmat[i][j] = 1. / (1. + r2 / D0_2) - 1/(1+4*Eff0);     // slightly extended to include those 8~9A
            //else Rmat[i][j] = -poor_penalty;
            //else Rmat[i][j] = -0.3;
        } else if(score_type == iGDT) {
            Rmat[i][j] = 0.;
            if(r2 > GDT_r2[3]){
                if(r2 < GDT_r2[3]*Eff0) Rmat[i][j] = 1.e-3;
                continue;
            }
            for(int m=3; m>=0; m--){
                if(r2 > GDT_r2[m]) break;
                Rmat[i][j] = 1. - m/4.;
            }
        } else {
            Rmat[i][j] = 1. / (1. + r2 / D0_2);
        }
    }
    }
}

void Salign::score_segment(int i, int j){
    auto& startA = pa->segstart;
    auto& midA = pa->segmid;
    auto& endA = pa->segend;
    auto& lenA = pa->seglen;

    auto& startB = pb->segstart;
    auto& midB = pb->segmid;
    auto& endB = pb->segend;
    auto& lenB = pb->seglen;

    double r2 = 0.; double rmsd=0.; double l1, d1, d2, d3, mind;
    int atom1_A_i, atom1_B_j, atom2_A_i, atom2_B_j, atom3_A_i, atom3_B_j;
    double my_d02 = 64.; //was 64 for non-ideal segs - early iteration may benefit from more flat score
    //double my_d02 = 81.; //was 64 for non-ideal segs - early iteration may benefit from more flat score
    double z_shift = 0.2; //0.3 was better but 0.2 consistent with SP-score

    atom2_A_i = midA[i];
    atom2_B_j = midB[j];
    d2 = xn1[atom2_A_i].distance2((*xbp_ideal)[atom2_B_j]);
    if (d2>100){Rmat_fragment[i][j]=-10.; return;} //TODO: optimize threshold for this shortcut heuristic
    d2 = (1./(1+d2/my_d02))-z_shift;

    atom1_A_i = startA[i];
    atom1_B_j = startB[j];

    d1 = xn1[atom1_A_i].distance2((*xbp_ideal)[atom1_B_j]);
    d1 = (1./(1+d1/my_d02))-z_shift;

    atom3_A_i = endA[i];
    atom3_B_j = endB[j];
    d3 = xn1[atom3_A_i].distance2((*xbp_ideal)[atom3_B_j]);
    d3 = (1./(1+d3/my_d02))-z_shift;

    Rmat_fragment[i][j] = (d1 + d2 + d3)/3.;
    l1 = min(lenA[i], lenB[j]); 
    Rmat_fragment[i][j] = max((Rmat_fragment[i][j])*pow(l1,segalpha), -10.); //d1w07a2, d1ivha1
}

void Salign::calRmat_fragment_bounded(int anchor[][3]){
    // NOTE RMAT_FRAGMENT IS NEVER EXPLICITY ZEROED OUT - MAY CAUSE PROBLEMS 
    auto& startA = pa->segstart;
    auto& midA = pa->segmid;
    auto& endA = pa->segend;
    auto& lenA = pa->seglen;

    auto& startB = pb->segstart;
    auto& midB = pb->segmid;
    auto& endB = pb->segend;
    auto& lenB = pb->seglen;

    double r2 = 0.; double rmsd=0.; double l1, d1, d2, d3, mind;
    int atom1_A_i, atom1_B_j, atom2_A_i, atom2_B_j, atom3_A_i, atom3_B_j;
    double my_d02 = 64.; //was 64 for non-ideal segs - early iteration may benefit from more flat score
    double z_shift = 0.2; //0.3 was better but 0.2 consistent with SP-score
    //Rmat_fragment.reset();
    Rmat_fragment.fill(-10.); //TODO: replace fill with smart method to avoid cases immediately overwritten
    for(int i=0; i<anchor[0][0]; i++)
    for(int j=0; j<anchor[1][0]; j++){
        if (pa->ssec[startA[i]]!=pb->ssec[startB[j]]){Rmat_fragment[i][j]=-100.; continue;}
        score_segment(i, j);
    }

    for (int z=0; z<3; z++){
        int i = anchor[0][z];
        int j = anchor[1][z];
        if ((i==-1) || (j==-1)) continue;
        score_segment(i,j);
    }

    for(int i=anchor[0][2]+1; i<nsA; i++)
    for(int j=anchor[1][2]+1; j<nsB; j++){
        if (pa->ssec[startA[i]]!=pb->ssec[startB[j]]){Rmat_fragment[i][j]=-100.; continue;}
        score_segment(i, j);
    }
}

void Salign::calRmat_ss(){
    auto& midA = pa->segmid;
    auto& midB = pb->segmid;

    for(int i=0; i<nsA; i++)
    for(int j=0; j<nsB; j++){
        if (pa->ssec[midA[i]]!=pb->ssec[midB[j]]){Rmat_fragment[i][j]=-10.; continue;}
        Rmat_fragment[i][j] = 1.;
    }
}


double Salign::optimize_score(int i0=0, int ie=6){
    double dmax1 = 0., d1;
    double r2 = powi(D0+2., 2);
//  int types[] = {0,2,2,1}; int rs[] = {0,r2+1,r2,0};
    int types[] = {0,1,2,2,1}; int rs[] = {0,0,r2+1,r2,0};
    for(int m=i0; m<ie; m++){
        int ierr = 0;
        //if(m >= 4) ierr = redo_fit(1,0); //20240124
        if(m >= 5) ierr = redo_fit(1,0); //20240124
        else ierr = redo_fit(types[m], rs[m]);
        if(ierr > 0 && m == 1){         // rare event in SP/GDT
            double dt = r2;
            for(int k=0; k<20; k++){
                dt *= 1.5;
                ierr = redo_fit(2, dt);
                if(ierr == 0) break;
            }
            if(ierr !=0) return dmax1;
        }
        if(ierr > 0) continue;
        if(ie>1 && m==0) continue;      // 2nd is always better
        d1 = calRscore();
        if(m>=4 && d1<dmax1) break; // mature break after 4 turns
        if(m>=4 && (d1-dmax1 < 0.05*dmax1)) break;  // mature break after 4 turns //5% convergence
        if(dmax1 < d1) dmax1 = d1;
    }
    return dmax1;
}

double Salign::optimize_align(){
    gap0 = gap1 = 0.;
    double dmax1=-100.;
    double gap0_dim[] = {0.5, final_gap0};
    bool broken = false;
    double d1;
    for(int k=0; k<2; k++){
        rotmol(nA, u_sv, (*xap), xn1);
        if (broken) break;
        gap0 = gap0_dim[k];
        for(int i=0; i<9; i++){
            if ((i==0) && (k==0)){
                calRmat();
                Align_DP();
                d1 = optimize_score(1,20);
            } else {
                init_block_align();
                refine_align(ialign);
                d1 = optimize_score(1,20);
            }
            if(d1-dmax1 < convergence_criterion*dmax1){
                break;
            }
            if(dmax1 < d1) dmax1 = d1;
        }
        if(dmax1 < d1) dmax1 = d1;
    }
    return dmax1;
}

void Salign::block_align(vector<int> ialign_cache[3]){
    //int i_end = 0;
    //int j_end = 0;
    int i_end = max(0, pa->segstart[ialign_cache[0][0]]-10);
    int j_end = max(0, pb->segstart[ialign_cache[1][0]]-10);
    int i_start, j_start;

    for (int z=0; z<ialign_cache[0].size(); z++){
        if ((ialign_cache[0][z]<0) || (ialign_cache[1][z]<0)) continue;
        i_start = pa->segstart[ialign_cache[0][z]];
        j_start = pb->segstart[ialign_cache[1][z]];

        //leading unaligned region trailing prior segment pair
        if (((i_start-i_end)!=0) && ((j_start-j_end)!=0)){ //may need to add gap penalty if only 1 has no connecting loop
            calRmat_bounded(i_end, i_start, j_end, j_start);
            Align_DP_bounded(i_end, i_start, j_end, j_start); 
        }
        i_end = pa->segend[ialign_cache[0][z]];
        j_end = pb->segend[ialign_cache[1][z]];
        
        //current segment pair
        calRmat_bounded(i_start, i_end, j_start, j_end); //for segments it might be better to use my sheet method - ie find best for first residue and then ungapped
        Align_DP_bounded(i_start, i_end, j_start, j_end);   
    }
    //i_start = pa->nres-1;
    //j_start = pb->nres-1;
    i_start = min(pa->nres-1, pa->segend[ialign_cache[0][ialign_cache[0].size()-1]]+10);
    j_start = min(pb->nres-1, pb->segend[ialign_cache[1][ialign_cache[0].size()-1]]+10);

    //trailing unaligned region
    if (((i_start-i_end)!=0) && ((j_start-j_end)!=0)){ 
        calRmat_bounded(i_end, i_start, j_end, j_start);
        Align_DP_bounded(i_end, i_start, j_end, j_start); 
    }
    traceback_DP_bounded();
}

void Salign::refine_align(vector<int> ialign_cache[3]){
    int i_end = 0;
    int j_end = 0;
    int i_start, j_start, i_match, j_match;

    //TODO: add leading tail
    int i_prev = -3;
    int j_prev = -3;

    for (int z=0; z<ialign_cache[0].size(); z++){
        if ((ialign_cache[0][z]<0) || (ialign_cache[1][z]<0)) continue;
        i_match = ialign_cache[0][z];
        j_match = ialign_cache[1][z];

        if ((i_match-i_prev)>2 && (j_match-j_prev)>2){
            
            if ((i_match-i_prev)>3 && (j_match-j_prev)>3){
                //leading loop
                //block 1
                i_start = max(i_prev+2+1, 0);
                i_end = min(i_match-1, nA-1);
                j_start = max(j_prev+1, 0);
                j_end = min(j_match-1, nB-1);

                calRmat_bounded(i_start, i_end, j_start, j_end);
                Align_DP_bounded(i_start, i_end, j_start, j_end);

                //block 2
                i_start = max(i_prev+1, 0);
                i_end = min(i_match-1, nA-1);
                j_start = max(j_prev+2+1, 0);
                j_end = min(j_match-1, nB-1);

                calRmat_bounded(i_start, i_end, j_start, j_end);
                Align_DP_bounded(i_start, i_end, j_start, j_end);

                //block 3 - must be last since it depends on block1 and 2 in align_dp
                i_start = max(i_prev+2+1, 0);
                i_end = min(i_match-1, nA-1);
                j_start = max(j_prev+2+1,0);
                j_end = min(j_match-1, nB-1);

                calRmat_bounded(i_start, i_end, j_start, j_end);
                Align_DP_bounded(i_start, i_end, j_start, j_end);           
            }

            // next match refinement
            //block 1
            i_start = max(i_match, 0);
            i_end = min(i_match+2, nA-1);
            j_start = max(j_match-2, 0);
            j_end = min(j_match-1, nB-1);

            calRmat_bounded(i_start, i_end, j_start, j_end);
            Align_DP_bounded(i_start, i_end, j_start, j_end);

            //block 2
            i_start = max(i_match-2, 0);
            i_end = min(i_match-1, nA-1);
            j_start = max(j_match, 0);
            j_end = min(j_match+2, nB-1);

            calRmat_bounded(i_start, i_end, j_start, j_end);
            Align_DP_bounded(i_start, i_end, j_start, j_end);

            //block 3 - must be last since it depends on block1 and 2 in align_dp
            i_start = max(i_match, 0);
            i_end = min(i_match+2, nA-1);
            j_start = max(j_match, 0);
            j_end = min(j_match+2, nB-1);

            calRmat_bounded(i_start, i_end, j_start, j_end);
            Align_DP_bounded(i_start, i_end, j_start, j_end);   

        } else {
            // next match refinement
            //block 1
            i_start = max(i_match+(3-(i_match-i_prev)), 0); //can probably merge with above
            i_end = min(i_match+2, nA-1);
            j_start = max(j_match-2, 0);
            j_end = min(j_match+1, nB-1);

            calRmat_bounded(i_start, i_end, j_start, j_end);
            Align_DP_bounded(i_start, i_end, j_start, j_end);

            //block 2
            i_start = max(i_match-2, 0);
            i_end = min(i_match+1, nA-1);
            j_start = max(j_match+(3-(j_match-j_prev)), 0);
            j_end = min(j_match+2, nB-1);

            calRmat_bounded(i_start, i_end, j_start, j_end);
            Align_DP_bounded(i_start, i_end, j_start, j_end);

            //block 3 - must be last since it depends on block1 and 2 in align_dp
            i_start = max(i_match+2, 0);
            i_end = min(i_match+2, nA-1);
            j_start = max(j_match+2, 0);
            j_end = min(j_match+2, nB-1);

            calRmat_bounded(i_start, i_end, j_start, j_end);
            Align_DP_bounded(i_start, i_end, j_start, j_end);               
        }
        i_prev = i_match;
        j_prev = j_match;
    }
    i_end = nA-1;
    j_end = nB-1;

    i_prev = max(i_prev, 0);
    j_prev = max(j_prev, 0);
    //rare case that seed alignment is corrupted by a bad segment match which leads to no residue score > 0

    //trailing unaligned region
    if (((i_end-i_prev)!=0) && ((j_end-j_prev)!=0)){ 
        calRmat_bounded(i_prev, i_end, j_prev, j_end);
        Align_DP_bounded(i_prev, i_end, j_prev, j_end); 
    }

    traceback_DP_bounded(); 
}

double Salign::optimize_block_align(){
    gap0 = gap1 = 0.;
    double dmax1=-100.;
    double gap0_dim[] = {0.5, final_gap0};
    vector<int> ialign_cache[3];
    for(int m=0; m<2; m++) ialign_cache[m] = ialign_fragment_sv[m];

    bool broken = false;
    double d1;
    
    for(int k=0; k<2; k++){
        rotmol(nA, u_sv, (*xap), xn1);
        if (broken) break;
        gap0 = gap0_dim[k];     
        for(int i=0; i<9; i++){
            if ((i==0) && (k==0)){
                init_block_align(); //may not be required first time
                block_align(ialign_cache);
                d1 = optimize_score(1,20); 
            } else {
                init_block_align();
                refine_align(ialign);
                d1 = optimize_score(1,20);
            }
            if(d1-dmax1 < convergence_criterion*dmax1){
                break;
            }
            if(dmax1 < d1) dmax1 = d1;
        }
        if(dmax1 < d1) dmax1 = d1;
    }
    return dmax1;
}



double Salign::run_Align1(){
    int ierr = 0;
    gap0 = gap1 = 0.;
    calRmat();
    Align_DP();
    return optimize_score(1);
}

double Salign::run_Align1_fragment(int anchor[][3]){ 
    int ierr = 0;
    double frag_score = 0.;

    //extend seed alignment to full structure
    if (pa->nseg<4 || pb->nseg<4){
        // This should never run
        rotmol(nA, u, (*xap), xn1); //rotate full set of coords 
        init_block_align();
        block_align(ialign_fragment);
        redo_fit(0);
        for(int m=0; m<riters; m++) redo_fit(1,1.);
        calRmat();
        frag_score = Align_DP();
    } else {
        //these are not relevant after scaling by lengths
        gap0 = 0.2;
        gap1 = 0.05;

        //may need to consider init_block_align - to reset seg3 alignment - maybe not
        calRmat_fragment_bounded(anchor);
        Align_DP_fragment_bounded(anchor);
        redo_fit_ideal(0);
        for(int m=0; m<riters; m++) redo_fit_ideal(1,1.); //this one might need even less

        calRmat_fragment_bounded(anchor);
        //second alignment may not be neccessary? re-use first and just score
        //can at least do local refinement rather than NxN

        frag_score = Align_DP_fragment_bounded(anchor);
    }

    Rscore = frag_score;
    return calRscore_frag();
}

void Salign::Run_all(){
    if(bscoreOnly >= 1){
        scoreOnly(); return;
    }
    Rscore_max = -100.;
    Rscore_sv = -100.;
    double d1;

    // Prefilter
    if ((pa->segstart.size()<1) || ((pb->segstart.size()<1))){
        double d1=1.;
    } else {
        calRmat_ss();
        gap0 = pref_gap0;
        gap1 = pref_gap1;
        if (bsingledom){
            d1 = Align_DP_fragment(true)/max(pa->nseg, pb->nseg);
        } else {
            d1 = Align_DP_fragment(true)/min(max(pa->nseg, pb->nseg),6);
        }
    }
    ss_prefilter = d1;

    if (ss_prefilter<ssprefcut){
        filtered = true;
        return;
    }

    // Main alignment
    if(nA>ma || nB>mb) alloc_large();
    // Not enough segments
    if ((pa->segstart.size()<4) || ((pb->segstart.size()<4))){
        Run3(); //fragment seeds (all atom)
        coarse = Rscore_max;
        if (Rscore_max>0){
            Rscore_sv = Rscore_max = -100; //my temporary hack due to difference in segment/all_atom score
            optimize_align();
        } else {
            filtered = true;
        }
    } else {
        RunSeg3(); //3segs
        RunSeg(); //2segs
        coarse = Rscore_max;
        if (coarse < coarsecut) { filtered=true; return; }
        if (Rscore_max>0){ //could easily have a more stringent filter here to greatly improve throughput
            Rscore_sv = Rscore_max = -100; //my temporary hack due to difference in segment/all_atom score
            optimize_block_align();
        } else {
            filtered = true;
        }
    }
    return;
}

double Salign::PrepareSeed(int i, int j, int offset_i, int offset_j, double rms0){
    auto& startA = pa->segstart;
    auto& midA = pa->segmid;
    auto& endA = pa->segend;
    auto& seglenA = pa->seglen;

    auto& startB = pb->segstart;
    auto& midB = pb->segmid;
    auto& endB = pb->segend;
    auto& seglenB = pb->seglen;

    for(int m=0; m<2; m++) ialign[m].clear();
    ialign[0].push_back(startA[i]);
    ialign[0].push_back(midA[i]);
    ialign[0].push_back(endA[i]);
    ialign[0].push_back(startA[i+offset_i]);
    ialign[0].push_back(midA[i+offset_i]);
    ialign[0].push_back(endA[i+offset_i]);
    ialign[1].push_back(startB[j]);
    ialign[1].push_back(midB[j]);
    ialign[1].push_back(endB[j]);
    ialign[1].push_back(startB[j+offset_j]);
    ialign[1].push_back(midB[j+offset_j]);
    ialign[1].push_back(endB[j+offset_j]);

    return redo_fit_ideal_seed(rms0);
}

void Salign::RunSeg(){
    auto& startA = pa->segstart;
    auto& midA = pa->segmid;
    auto& endA = pa->segend;
    auto& seglenA = pa->seglen;

    auto& startB = pb->segstart;
    auto& midB = pb->segmid;
    auto& endB = pb->segend;
    auto& seglenB = pb->seglen;

    if (startA.size()<2) return;
    if (startB.size()<2) return;
    double rms1;
    double rms0 = segcut;

    double opt_rms = 1000;
    int opt_i = 0;
    int opt_j = 0;
    int opt_shift_i = 1;
    int opt_shift_j = 1;
    int anchor[2][3];

    for(int i=0; i<startA.size()-1; i+=1){
        for(int j=0; j<startB.size()-1; j+=1){
        for (int skip_offset_i=1; skip_offset_i<3; skip_offset_i+=1){
        for (int skip_offset_j=1; skip_offset_j<3; skip_offset_j+=1){
            if (((i+skip_offset_i+1)>startA.size()) || ((j+skip_offset_j+1)>startB.size())) continue;
            if (pa->ssec[startA[i]]!=pb->ssec[startB[j]]) continue;
            if (pa->ssec[startA[i+skip_offset_i]]!=pb->ssec[startB[j+skip_offset_j]]) continue;
            if (success3.find(make_tuple(i, i+skip_offset_i, j, j+skip_offset_j)) != success3.end()) continue;
            seeds++;
            rms1 = PrepareSeed(i, j, skip_offset_i, skip_offset_j, rms0);
            if (rms1 < opt_rms){
                opt_rms = rms1;
                opt_i = i;
                opt_j = j;
                opt_shift_i = skip_offset_i;
                opt_shift_j = skip_offset_j;
            }
            if(rms1 > rms0) continue;
            valid_seeds++;
            for(int m=0; m<riters; m++) redo_fit_ideal(1,1.);
            
            anchor[0][0] = i;
            anchor[1][0] = j;
            if (skip_offset_i>1) {anchor[0][1] = i+1;}
            else {anchor[0][1] = -1;}
            if (skip_offset_j>1) anchor[1][1] = j+1;
            else {anchor[1][1] = -1;}
            anchor[0][2] = i+skip_offset_i;
            anchor[1][2] = j+skip_offset_j;
            double d1 = run_Align1_fragment(anchor);
        }
        }
        }
    }
    if (opt_rms < 8. && opt_rms>=rms0){
        seeds++;
        valid_seeds++;
        rms1 = PrepareSeed(opt_i, opt_j, opt_shift_i, opt_shift_j, 1000.);
        for(int m=0; m<riters; m++) redo_fit_ideal(1,1.);
        anchor[0][0] = opt_i;
        anchor[1][0] = opt_j;
        if (opt_shift_i>1) {anchor[0][1] = opt_i+1;}
        else {anchor[0][1] = -1;}
        if (opt_shift_j>1) anchor[1][1] = opt_j+1;
        else {anchor[1][1] = -1;}
        anchor[0][2] = opt_i+opt_shift_i;
        anchor[1][2] = opt_j+opt_shift_j;
        double d1 = run_Align1_fragment(anchor);
    } 
}

void Salign::RunSeg3(){
    auto& startA = pa->segstart;
    auto& midA = pa->segmid;
    auto& endA = pa->segend;
    auto& seglenA = pa->seglen;

    auto& startB = pb->segstart;
    auto& midB = pb->segmid;
    auto& endB = pb->segend;
    auto& seglenB = pb->seglen;

    if (startA.size()<3) return;
    if (startB.size()<3) return;
    double rms1;
    double rms0 = 6; //6 is good
    int anchor[2][3];
    for(int i=0; i<startA.size()-2; i+=1){
        for(int j=0; j<startB.size()-2; j+=1){
            if (pa->ssec[startA[i]]!=pb->ssec[startB[j]]) continue;
            if (pa->ssec[startA[i+1]]!=pb->ssec[startB[j+1]]) continue;
            if (pa->ssec[startA[i+2]]!=pb->ssec[startB[j+2]]) continue;
            seeds++;
        
            for(int m=0; m<2; m++) ialign[m].clear();
            
            ialign[0].push_back(startA[i]);
            ialign[0].push_back(midA[i]);
            ialign[0].push_back(endA[i]);
            ialign[0].push_back(startA[i+1]);
            ialign[0].push_back(midA[i+1]);
            ialign[0].push_back(endA[i+1]);
            ialign[0].push_back(startA[i+2]);
            ialign[0].push_back(midA[i+2]);
            ialign[0].push_back(endA[i+2]);
            ialign[1].push_back(startB[j]);
            ialign[1].push_back(midB[j]);
            ialign[1].push_back(endB[j]);
            ialign[1].push_back(startB[j+1]);
            ialign[1].push_back(midB[j+1]);
            ialign[1].push_back(endB[j+1]);
            ialign[1].push_back(startB[j+2]);
            ialign[1].push_back(midB[j+2]);
            ialign[1].push_back(endB[j+2]);

            rms1 = redo_fit_ideal_seed(rms0);
            if(rms1 > rms0) continue;
            valid_seeds++;
            for(int m=0; m<riters; m++) redo_fit_ideal(1,1.); //this wasnt 1 - check
            
            anchor[0][0] = i;
            anchor[0][1] = i+1;
            anchor[0][2] = i+2;
            anchor[1][0] = j;
            anchor[1][1] = j+1;
            anchor[1][2] = j+2;

            double d1 = run_Align1_fragment(anchor);
            success3.insert(make_tuple(i, i+1, j, j+1));
            success3.insert(make_tuple(i+1, i+2, j+1, j+2));
            //gapped
            success3.insert(make_tuple(i, i+2, j, j+2));
        }
    }
}


// init seed by fragment
void Salign::Run3(){
    Rscore_sv = -100.;
    int nmin = min(nA, nB);
    int nfrag = max(20, min(50, nmin/5));
    if(fragsize >= 1) nfrag = max(fragsize, 5);
    if(nfrag > nmin) nfrag = nmin;
    double rms1;
    int n0 = nfrag, n1 = nfrag;
    double rms0 = pow(nfrag, 1./3) * 1.2;
    int opt_i, opt_j;
    double opt_rms = 1000;

    for(int i=0; i<nA-nfrag+1; i++)
    for(int j=0; j<nB-nfrag+1; j+=n1){
        for(int m=0; m<2; m++) ialign[m].clear();
        for(int m=0; m<nfrag; m++){
            ialign[0].push_back(i+m);
            ialign[1].push_back(j+m);
        }
        redo_fit(0);
        calrms(rms1);
        if(rms1 < opt_rms){
            opt_rms = rms1;
            opt_i = i;
            opt_j = j;
        }
        if(rms1 > rms0) continue;
        for(int m=0; m<3; m++) redo_fit(1, 1.);
        double d1 = run_Align1();
    }
    if (opt_rms>rms0){
        for(int m=0; m<2; m++) ialign[m].clear();
        for(int m=0; m<nfrag; m++){
            ialign[0].push_back(opt_i+m);
            ialign[1].push_back(opt_j+m);
        }
        redo_fit(0);
        calrms(rms1);
        for(int m=0; m<3; m++) redo_fit(1, 1.);
        double d1 = run_Align1();
    }
}

// calculate score based on values of xn1/xbp; ialign & D0_2
double Salign::calRscore(){
    Rscore = 0.;
    for(int i=0; i<ialign[0].size(); i++){
        int i1 = ialign[0][i], i2 = ialign[1][i];
        if(i1 < 0 || i2 < 0) continue;
        double r2 = xn1[i1].distance2((*xbp)[i2]);
        if(score_type == iSP) {
            if(r2 > D0_2*4) continue;
            Rscore += 1.25 * (1./ (1. + r2/D0_2) - 0.2);
        } else if(score_type == iGDT) {
            double dt = 0.;
            for(int m=3; m>=0; m--){
                if(r2 > GDT_r2[m]) break;
                dt = 1 - m/4.;
            }
            Rscore += dt;
        } else {
            Rscore += 1. / (1. + r2/D0_2);
        }
    }
    save_align();
    return Rscore;
}

double Salign::calRscore_frag(){
    save_align();
    return Rscore;
}

void Salign::calscores_all(){
    double nali, rms, GDT, SP0, LG0, nid1;
    SP0 = LG0 = GDT = nali = rms = nid1 = 0.;
    double D2 = D00*D00;
// TMs
    double TMs[3], lTM[3], DTM2[3];
    bzero(TMs, sizeof(TMs));
    lTM[0] = (nA + nB) * 0.5;
    lTM[1] = min(nA, nB);
    lTM[2] = max(nA, nB);
    for(int m=0; m<3; m++) {
        double dt = getd0_TM(lTM[m]);
        DTM2[m] =  dt*dt;
    }
//
    for(int i=0; i<ialign[0].size(); i++){
        int i1 = ialign[0][i], i2 = ialign[1][i];
        if(i1 < 0 || i2 < 0) continue;
        double r2 = xn1[i1].distance2((*xbp)[i2]);
        for(int m=0; m<3; m++) TMs[m] += 1. / (1. + r2/DTM2[m]);
        LG0 += 1. / (1. + r2/D2);
        for(int m=0; m<4; m++){
            if(r2 > GDT_r2[m]) continue;
            GDT += 1 - m/4.; break;
        }
        if(r2 < D2*4) SP0 += 1.25 * (1/ (1. + r2/D2) - 0.2);
        if(r2 < 64.) { nali ++; rms += r2; }
        if(pa->resid[i1] == pb->resid[i2]) nid1 ++;
    }
    scores_all[ieLA] = nali;
    nali = max(1., nali); rms = sqrt(rms / nali);
    scores_all[ieRMS] = rms;
    scores_all[ieSP] = SP0 / Lmin;
    if(score_type == iSP){
        int num[2][4]; calLali2(&num[0][0], Denv);
        assert(num[0][0] == num[1][0]);
        scores_all[ieLE] = (num[0][1] + num[1][1])*0.5 + num[0][0];
    }
    for(int m=0; m<3; m++){
        TMs[m] /= lTM[m];
    }
    scores_all[ieTMa] = TMs[0];
    scores_all[ieTMb] = TMs[1];
    scores_all[ieTMc] = TMs[2];
    scores_all[ieGDT] = GDT / Lmin;
    scores_all[ieLG] = LG0 / Lmin;
    scores_all[ieSEQ] = 100.*nid1 / Lmin;
    return;
}
bool Salign::save_align(){
    if(Rscore_sv > Rscore) return -1;
    Rscore_sv = Rscore;
    for(int m=0; m<12; m++) u_sv[0][m] = u[0][m];
    for(int m=0; m<2; m++) ialign_sv[m] = ialign[m];
    for(int m=0; m<2; m++) ialign_fragment_sv[m] = ialign_fragment[m];
//
    if(Rscore_max > Rscore) return 0;
    Rscore_max = Rscore;
    for(int m=0; m<12; m++) u_max[0][m] = u[0][m];
    for(int m=0; m<2; m++) ialign_max[m] = ialign[m];
    for(int m=0; m<2; m++) ialign_fragment_max[m] = ialign_fragment[m];
    return 1;
}
void Salign::print_max(FILE *fp){
    string sinfo;
    print_max(sinfo);
    fprintf(fp, "%s", sinfo.c_str());
}

void Salign::restore_max(){
    for(int m=0; m<12; m++) u[0][m] = u_max[0][m];
    for(int m=0; m<2; m++) ialign[m] = ialign_max[m];
};

void Salign::print_max(string &sinfo){
    if (Rscore_max>0){
        for(int m=0; m<12; m++) u[0][m] = u_max[0][m];
        for(int m=0; m<2; m++) ialign[m] = ialign_max[m];
        score_final();
    }
    prtali(sinfo);
}
void Salign::pro_swap(int bmat=0){
    breverse = (! breverse);
    Protein2 *pt = pa; pa = pb; pb = pt;
    nA = pa->nres; nB = pb->nres;
    xap = &pa->x; xbp = &pb->x;
    if(! bmat) return;
// swap the matrix & ialign
    for(int i=0; i<3; i++)
    for(int j=i+1; j<3; j++){
        double dt = u[i][j];
        u[i][j] = u[j][i]; u[j][i] = dt;
    }
    double dat[3];
    for(int m=0; m<3; m++) dat[m] = -dot_product(u[m], u[3]);
    for(int m=0; m<3; m++) u[3][m] = dat[m];
    vector<int> vt = ialign[0];
    ialign[0] = ialign[1]; ialign[1] = vt;
}

double Salign::run_pairwise(Protein2 *a, Protein2 *b){
    init(a, b);
    Run_all();
    if (getRmax()<0) 
        return 0.;
    restore_max();
    score_final();
    return calSPscore();
}

int Salign::init(Protein2 *a, Protein2 *b, int brev){       // nA < nB
    Rscore_max = Rscore_sv = -100;
    pa = a; pb = b; breverse = 0;
    if(brev && pa->nres > pb->nres){
        pb = a; pa = b; breverse = 1;
    }
    nA = pa->nres; nB = pb->nres;
    nsA = pa->nseg; nsB = pb->nseg;
    xap = &pa->x; xbp = &pb->x;
    xap_ideal = &pa->x_ideal; xbp_ideal = &pb->x_ideal;
    int nmax = max(nA, nB);
    list2[0].reserve(nmax); //might be useless
    list2[1].reserve(nmax);
    wfit.reserve(nmax);
    u_sv[0][0] = u_max[0][0] = 1.;
    u_sv[1][1] = u_max[1][1] = 1.;
    u_sv[2][2] = u_max[2][2] = 1.;
    if(nA < 5){
        //fprintf(stderr, "empty file: %s\n", pa->getname().c_str());
        return 1;
    }
    if(nB < 5) {
        //fprintf(stderr, "empty file: %s\n", pb->getname().c_str());
        return 1;
    }
    Lmin = min(nA, nB);
    D0 = D00;
    if(score_type == iTM) D0 = getd0_TM(double(Lmin));
    D0_2 = D0 * D0;
//
    if(nsA>msa || nsB>msb) alloc_x();
    if(nA>ma || nB>mb) alloc_lengths();
    if(DEBUG > 0) printf("Nres: %s %d -- %s %d; cutoff: %.1f\n", 
        pa->name.c_str(), nA, pb->name.c_str(), nB, D0);
    return 0;
}
// swap back if in reverse status, and compute all scores for Rscore_max
void Salign::score_final(){
    if(breverse) pro_swap(1);
    rotmol(nA, u, (*xap), xn1);
    calscores_all();
}
void Salign::calLali2(int *num, double R0){
    vector<int> idx[2];
    idx[0].resize(nA); idx[1].resize(nB);
    for(int i=0; i<nA; i++) idx[0][i] = -1;
    for(int i=0; i<nB; i++) idx[1][i] = -1;
//
    vector<int> aligned1;
    vector<int> aligned2;
    for(int i=0; i<ialign[0].size(); i++){
        int i1 = ialign[0][i], i2 = ialign[1][i];
        if(i1 < 0 || i2 < 0) continue;
        double r2 = xn1[i1].distance2((*xbp)[i2]);
        if(r2 > 64.) continue;
        idx[0][i1] = idx[1][i2] = 1;        // core region
        aligned1.push_back(i1);
        aligned2.push_back(i2);
    }
    pa -> calnneib(idx[0], aligned1, R0);
    pb -> calnneib(idx[1], aligned2, R0);
    for(int m=0; m<8; m++) num[m] = 0;
    for(int k=0; k<2; k++){
        for(int i=0; i<idx[k].size(); i++){
            int it = idx[k][i] - 1;
            if(it < 0 || it > 3) continue;
            num[k*4 + it] ++;
        }
    }
}
void Salign::delete_x(){
    if(ma<1 && mb<1) return;
}
void Salign::alloc_x(){
    //delete_x();
    msa = nsA; msb = nsB;
    Idir_fragment.resize(nsA+1, nsB+1);
    Rmat_fragment.resize(nsA+1, nsB+1);
    Smat_fragment.resize(nsA+1, nsB+1);
}

void Salign::alloc_lengths(){
    //delete_x();
    int nmax = max(nA, nB);
    wfit.resize(nmax, 1.);
    xn1.resize(nmax); xn2.resize(nmax);
}

void Salign::alloc_large(){
    //delete_x();
    ma = nA; mb = nB;
    Idir.resize(ma+1, mb+1);
    Rmat.resize(ma+1, mb+1);
    Smat.resize(ma+1, mb+1);
}
void Salign::prtali(FILE *fp){
    string sinfo;
    prtali(sinfo);
    fprintf(fp, "%s", sinfo.c_str());
}
void Salign::calSPscores(double *sp){
    double sp0 = scores_all[ieSP]*Lmin;
    double L1 = max(scores_all[ieLE], 3.); //Note: fractional length possible
    //sp[0] = sp0 / (pow(L1, 1.-Alpha) * 3.75);
    sp[0] = sp0 / (pow(L1, 1.-Alpha) * scaling_factor);
    //sp[1] = sp0 / (pow(nA, 1.-Alpha) * 3.75);
    //sp[1] = sp0 / max(nA,nB);
    sp[1] = sp0;
//  sp[2] = sp0 / (pow(nB, 1.-Alpha) * 3.75);
//  sp[1] = sp0 / (pow((nA+nB)*0.5, 1.-Alpha) * 3.75);
//  sp[2] = sp0 / (pow(Lmin, 1.-Alpha) * 3.75);
    //sp[3] = 100. / (1. + exp(-(sp[0] - 0.523) / 0.044));
    sp[2] = ss_prefilter;
    //printf("SP: %.3f\n", sp[0]);
}
void Salign::prtdis_pair(string &sinfo){
    sinfo += "#distance pair\n";
    char str[201];
    for(int i=0; i<ialign[0].size(); i++){
        int i1 = ialign[0][i], i2 = ialign[1][i];
        if(i1 < 0 || i2 < 0) continue;
        double r2 = xn1[i1].distance2((*xbp)[i2]);
        char c1=rnam1_std[pa->resid[i1]], c2 = rnam1_std[pb->resid[i2]];
        sprintf(str, "%d %d %c%c %.1f\n", pa->seq0[i1], pb->seq0[i2], c1, c2, sqrt(r2));
        sinfo += str;
    }
}
void Salign::prtali(string &sinfo){
    sinfo = "";
    char str[max(nA+nB+100,2001)];
//
    double SPs[4];      // SPe, SPa, SPb
    if(score_type == iSP) {
        calSPscores(SPs);
        if (SPs[0]<reportcutoff) return;
    }
//
    if(iprint >= 99) {
        sprintf(str, "################## Alignment Report ##################\n"); sinfo += str;
        sprintf(str, "%s %s: Length= %d %d\n", pa->name.c_str(), pb->name.c_str(), nA, nB); sinfo += str;
        if(score_type == iSP){
            double p0 = SPs[3];
            sprintf(str, "Pfold= %.1f %%; SPe/SPa/SPb= %.3f %.3f %.3f ;Effective_Length: %d\n", p0, SPs[0], SPs[1], SPs[2], int(scores_all[ieLE])); sinfo += str;
        }
        sprintf(str, "RMSD/Nali= %.2f / %d ;GDT= %.3f ;TMscore(a,b,c)= %.3f %.3f %.3f SEQID= %.1f%%\n", scores_all[ieRMS], int(scores_all[ieLA]+0.1), scores_all[ieGDT], scores_all[ieTMa], scores_all[ieTMb], scores_all[ieTMc], scores_all[ieSEQ]); sinfo += str;
        sprintf(str, "\nRotation Matrix:\n"); sinfo += str;

    } else {
        double e1 = scores_all[score_type];
        if(score_type == iSP) e1 = SPs[0];
        if(e1 < cutoff) return;
//
        sprintf(str, "%s %s ", pa->name.c_str(), pb->name.c_str()); sinfo+=str;
        if(score_type == iSP){
            sprintf(str, "%.3f %.3f %.3f %d %d %d %.1f ", SPs[0], SPs[1], SPs[2], nA, nB, int(scores_all[ieLE]), scores_all[ieSEQ]);
        } else if(score_type == iTM){
            sprintf(str, "%.3f %.3f %.3f ", scores_all[ieTMa], scores_all[ieTMb], scores_all[ieTMc]);
        } else {
            sprintf(str, "%.3f ", scores_all[score_type]);
        }
        sinfo += str;
        sprintf(str, "%d %d %d %.3f\n", int(scores_all[ieLA]+0.1), seeds, valid_seeds, coarse); sinfo+=str;
    }
    if(iprint <= 1) return;
//
// print matrix
    for(int m=0; m<3; m++){
        sprintf(str, "%10.5f %8.5f %8.5f %8.5f\n", u[3][m], u[m][0], u[m][1], u[m][2]); sinfo += str;
    }
    if(iprint == 2) return;
//
    string sa1, sa2, sdis, dis_info=""; sa1 = sa2 = sdis = "";
    for(int i=0; i<ialign[0].size(); i++){
        int i1=ialign[0][i], i2=ialign[1][i];
        if(i1 < 0) sa1 += '-';
        else sa1 += rnam1_std[pa->resid[i1]];
        if(i2 < 0) sa2 += '-';
        else sa2 += rnam1_std[pb->resid[i2]];
        char c1 = ' ';
        if(i1>=0 && i2>=0) {
            double r2 = xn1[i1].distance2( (*xbp)[i2] );
            if(r2 <= 16) c1 = ':';
            else if(r2 < 25) c1 = ';';
            else if(r2 < 64) c1 = '.';
        }
        sdis += c1;
        if(iprint > 999) {
            double r2 = xn1[i1].distance2( (*xbp)[i2] );
            sprintf(str, "%c%d %c%d %.2f\n", rnam1_std[pa->resid[i1]], i1, rnam1_std[pb->resid[i2]], i2, sqrt(r2));
            dis_info += str;
        }
    }
//
    int is1=0, is2=0, ie1=0, ie2=0;
    for(int i=0; i<ialign[0].size(); i++){
        if(ialign[0][i] >= 0) {is1=ialign[0][i]; break;}
    }
    for(int i=ialign[0].size()-1; i>=0; i--){
        if(ialign[0][i] >= 0) {ie1=ialign[0][i]; break;}
    }
    for(int i=0; i<ialign[0].size(); i++){
        if(ialign[1][i] >= 0) {is2=ialign[1][i]; break;}
    }
    for(int i=ialign[1].size()-1; i>=0; i--){
        if(ialign[1][i] >= 0) {ie2=ialign[1][i]; break;}
    }
    if(iprint >= 99){
        sprintf(str, "\nAlignment: %s %s\n", pa->name.c_str(), pb->name.c_str()); sinfo += str;
        sprintf(str, "(':' denotes the residue pairs of distance <= 4A, and '.' denotes <=8A)\n"); sinfo += str;
    }
    if(bfullalign){
        while(is1>0 || is2>0){
            char c1 = '-', c2 = '-';
            if(is1 > 0) c1 = rnam1_std[pa->resid[--is1]];
            if(is2 > 0) c2 = rnam1_std[pb->resid[--is2]];
            sa1 = c1 + sa1;
            sa2 = c2 + sa2;
            sdis = ' ' + sdis;
        }
        while(ie1<nA-1 || ie2<nB-1){
            char c1 = '-', c2 = '-';
            if(ie1 < nA-1) c1 = rnam1_std[pa->resid[++ie1]];
            if(ie2 < nB-1) c2 = rnam1_std[pb->resid[++ie2]];
            sa1 = sa1 + c1;
            sa2 = sa2 + c2;
            sdis = sdis + ' ';
        }
        sprintf(str, "%s\n", sa1.c_str()); sinfo += str;
        sprintf(str, "%s\n", sdis.c_str()); sinfo += str;
        sprintf(str, "%s\n", sa2.c_str()); sinfo += str;
    } else {
        sprintf(str, "%-4d %s %d\n", is1, sa1.c_str(), ie1); sinfo += str;
        sprintf(str, "%-4s %s\n", "", sdis.c_str()); sinfo += str;
        sprintf(str, "%-4d %s %d\n", is2, sa2.c_str(), ie2); sinfo += str;
    }
    sinfo += dis_info;
    if(iprint > 200) prtdis_pair(sinfo);
}
//
int Salign::fit1(vector<int> *iali2){
    double xc1[3], xc2[3];
    int nfit = 0;
    for(int i=0; i<iali2[0].size(); i++){
        int i1=iali2[0][i], i2=iali2[1][i];
        if(i1<0 || i2<0) continue;
        wfit[nfit] = 1.;
        for(int m=0; m<3; m++){
            xn1[nfit][m] = (*xap)[i1][m];
            xn2[nfit][m] = (*xbp)[i2][m];
        }
        nfit ++;
    }
    if(nfit < 3) die("too small alignment size: %d", nfit);
    double ds = 0.;
    for(int i=0; i<nfit; i++) ds += wfit[i];
    if(ds > 1.0e-8){
        for(int i=0; i<nfit; i++) wfit[i] /= ds;
    } else {
        for(int i=0; i<nfit; i++) wfit[i] = 1./nfit;
    }
    translate2(nfit, wfit, xn1, xn1, xc1);
    translate2(nfit, wfit, xn2, xn2, xc2);
    quatfit(nfit, wfit, xn2, xn1, u);
    for(int m=0; m<3; m++) u[3][m] = xc2[m] - dot_product(u[m], xc1);
    rotmol(nA, u, (*xap), xn1);
    return 0;
}
void Salign::scoreOnly(){
    Rscore_max = Rscore_sv = -100;
    for(int m=0; m<2; m++) ialign[m].clear();
    assert(nA <= nB);
    if(bscoreOnly == 1){
        map<string, int> resnum2id;
        for(int i=0; i<pa->resnum.size(); i++){
            resnum2id[pa->resnum[i]] = i;
        }
        for(int i=0; i<pb->resnum.size(); i++){
            if(resnum2id.count(pb->resnum[i]) <= 0) continue;
            ialign[0].push_back(resnum2id[pb->resnum[i]]);
            ialign[1].push_back(i);
        }
    } else {
        for(int i=0; i<nA; i++){
            ialign[0].push_back(i);
            ialign[1].push_back(i);
        }
    }
    int nfit0=ialign[0].size(),  nfit=nfit0;
    if(nfit0 < 3) die("too small number of match residues for score: %d", nfit0);
    fit1(ialign);
    double rms=0; calrms(rms);
    if(iprint>=99) printf("initial RMSD: %.2f / %d\n", rms, nfit);
//
    vector<int> iali2[2];
    for(int k=0; k<10; k++){
        for(int i=0; i<nfit0-nfit+1; i++){
            for(int m=0; m<2; m++) iali2[m].clear();
            for(int m=0; m<nfit; m++){
                iali2[0].push_back( ialign[0][i+m] );
                iali2[1].push_back( ialign[1][i+m] );
            }
            fit1(iali2);
            optimize_score(1);
        }
        if(nfit <= 4) break;
        nfit /= 4;
        if(nfit < 4) nfit = 4;
    }
    Rscore = Rscore_max;
}
