#include "protein2.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <stdint.h>

void Protein2::rdpdb(string fn){
    info_err = "";
    FILE *fp = openfile(fn, "r");
    if(fp == NULL) {fprintf(stderr, "wrong file: %s\n", fn.c_str()); return;}
    char str[121], rn[8], an[5]; string rinfo0="", rinfo;
    double xt[3];
    name = basename(fn.c_str());
    while(fgets(str,120,fp) != NULL){
        if(strstr(str, "END") == str) break;
        if(strstr(str, "TER") == str) break;
        if(strstr(str, "ATOM ") != str) continue;
        sscanf(str+13, "%3s", an); an[2] = '\0';
        if(strcmp(an, "CA") != 0) continue;
        sscanf(str+17, "%3s", rn); rn[3] = '\0';
        string line(str);
        rinfo = line.substr(17, 10);
        if(rinfo == rinfo0) continue;
        for(int m=0; m<3; m++) xt[m] = strtod(str+30+8*m, NULL);
        x.push_back(Xvec(xt));
        int id = aaDefine(rn, DEBUG);
        seq0.push_back(strtol(str+22,NULL,10));
        resid.push_back(id);
        sscanf(str+22, "%7s", rn); rn[7] = '\0';
        resnum.push_back(rn);
        rinfo0 = rinfo;
    }
    nres = resid.size();
    fclose(fp);
    if(nres < 3) {
        info_err += " few residues";
        return;
    }
    calSS();
}


void Protein2::rdideal(string fn){
    char line[400];
    char aa;
    int ss;
    double coord[3];
    double ideal_coord[3];
    int start, mid, end, tmp_nres, tmp_nseg;
    vector<int> ssec0;
    name = basename(fn.c_str());
    FILE* fp = fopen(fn.c_str(), "r");
    if(fp == NULL) {fprintf(stderr, "wrong file: %s\n", fn.c_str()); return;}

    fgets(line, 400, fp);
    sscanf(line, "%d %d", &nres, &nseg);


    for (int i=0; i<nres; i++){
        fgets(line, 400, fp);
        sscanf(line, "%c %lf %lf %lf %d %lf %lf %lf",
            &aa, &coord[0], &coord[1], &coord[2], &ss, &ideal_coord[0], &ideal_coord[1], &ideal_coord[2]);
        int id = aaDefine1(aa); 
        resid.push_back(id);
        x.push_back(Xvec(coord));
        x_ideal.push_back(Xvec(ideal_coord));
        ssec0.push_back(ss);
    }

    for (int i=0; i<nseg; i++){
        fgets(line, 400, fp);
        sscanf(line, "%d %d %d", &start, &mid, &end);
        segstart.push_back(start);
        segmid.push_back(mid);
        segend.push_back(end);
        seglen.push_back((end-start)+1);
    }

    fclose(fp);
    ssec = ssec0;
}


Protein2::Protein2(int nr, string seq1, vector<vector<double> > xp){
    info_err = "";
    name = "UNK";
    nres = nr;
    resid.resize(nr);
    for(int i=0; i<nr; i++) resid[i] = aaDefine1(seq1[i]);
    x.resize(nr);
    for(int i=0; i<nr; i++){
        x[i].setx(&xp[i][0]);
    }
    if(nres < 3) {
        info_err += " few residues";
        return;
    }
}
void Protein2::calnneib(vector<int> &idx, vector<int> &aligned, double R0){
    const double R02 = R0*R0;
    for(int i=0; i<nres; i++){
        if(idx[i] == 1) continue; //only consider non-core residues
        for (const int& j: aligned){
            const double r2 = distance2(i, j);
            if(r2 < R02){
                idx[i] = 2; 
                break;
            }
        }
    }
    return;
}
void Protein2::calnneib0(vector<int> &idx){
    assert(idx.size() == nres);
    double r2_dim[] = {10., 12., 14.};
    for(int i=0; i<nres; i++){
        if(idx[i] == 1) continue;
        double rmin2 = 10000.;
        for(int j=0; j<nres; j++){
            if(idx[j] != 1) continue;
            double r2 = distance2(i, j);
            if(rmin2 < r2) continue;
            rmin2 = r2;
            if(rmin2 <= r2_dim[0]*r2_dim[0]) break;
        }
        for(int m=0; m<3; m++){
            if(rmin2 > r2_dim[m]*r2_dim[m]) continue;
            idx[i] = 2 + m; break;
        }
    }
}
void Protein2::calSS2(vector<int> &ssec0){
    ssec = ssec0;
    nres = ssec.size();

    vector<int> ifrag[2];
    int i1 = -1;
    for(int i=0; i<nres; i++){
        if(ssec[i]>0 && ssec[i]!=i1){
            ifrag[0].push_back(i); 
        } 
        if (ssec[i]!=i1 && i1>0){
            ifrag[1].push_back(i-1); //after hitting a new type, add the previous residue to segment end
        }
        i1 = ssec[i];
    }
    if(ifrag[1].size() < ifrag[0].size()) ifrag[1].push_back(nres-1);

    for(int i=0; i<ifrag[0].size(); i++){
        int i0=ifrag[0][i], i1=ifrag[1][i];
        int len = 1+i1-i0;
        bool strand = ((len >= 3) && (ssec[i0]==2));
        bool helix = ((len >= 6) && (ssec[i0]==1)); //was 10 originally
        if((strand) || (helix)) {
            segstart.push_back(i0);
            segmid.push_back(i0+((len)/2));
            segend.push_back(i1);
            seglen.push_back(len);
        }
    }
    nseg = segstart.size();
}

void Protein2::add_ideal(int nr, vector<vector<double> > xp){
    x_ideal.resize(nr);
    for(int i=0; i<nr; i++){
        x_ideal[i].setx(&xp[i][0]);
    }
}

void Protein2::add_dummy_ideal(){
    x_ideal.resize(nres);
    for(int i=0; i<nres; i++){
        x_ideal[i].setx(&x[i][0]);
    }
}

void Protein2::calSS(){
    double Rdis[nres][3];
    const double Ra[] = {5.45, 5.18, 6.37}, Da = 2.1;
    const double Rb[] = {6.1, 10.4, 13}, Db = 1.42;
    if(nres <= 3) {
        fprintf(stderr, "%s: few residues %d", name.c_str(), nres);
        return;
    }
    ssec.resize(nres);
    for(int k=0; k<3; k++){
        for(int i=0; i<nres-k-2; i++) Rdis[i][k] = sqrt(distance2(i, i+k+2));
    }
    for(int i=0; i<nres; i++) ssec[i] = 0;
    for(int i=0; i<nres-5; i++){
        int iss = 0;
        if(fabs(Rdis[i][2]-Ra[2])<Da &&
            fabs(Rdis[i][1]-Ra[1])<Da && fabs(Rdis[i+1][1]-Ra[1])<Da &&
            fabs(Rdis[i][0]-Ra[0])<Da && fabs(Rdis[i+1][0]-Ra[0])<Da &&
            fabs(Rdis[i+2][0]-Ra[0])<Da) iss = 1;

        else if(fabs(Rdis[i][2]-Rb[2])<Db &&
            fabs(Rdis[i][1]-Rb[1])<Db && fabs(Rdis[i+1][1]-Rb[1])<Db &&
            fabs(Rdis[i][0]-Rb[0])<Db && fabs(Rdis[i+1][0]-Rb[0])<Db &&
            fabs(Rdis[i+2][0]-Rb[0])<Db) iss = 2;
        else continue;
        ssec[i+2] = iss;
    }

    ssec[0] = ssec[1] = ssec[2];
    ssec[nres-1] = ssec[nres-2] = ssec[nres-3];

// fragment
    vector<int> ifrag[2];
    int i1 = 0;
    if (ssec[0]>0) {ifrag[0].push_back(0); i1=ssec[0];}
    for(int i=1; i<nres; i++){
        if(ssec[i]>0 && ssec[i-1]!=ssec[i]){
            ifrag[0].push_back(i); i1 = ssec[i];
        } else if(ssec[i] != i1){
            if(ifrag[1].size() < ifrag[0].size()) ifrag[1].push_back(i);
            i1 = ssec[i];
        }
    }
    if(ifrag[1].size() < ifrag[0].size()) ifrag[1].push_back(nres);

// head & tail
    for(int i=0; i<ifrag[0].size(); i++){
        int i0=ifrag[0][i], i1=ifrag[1][i];
        bool strand = ((i1-i0 >= 3) && (ssec[i0]==2));
        bool helix = ((i1-i0 >= 6) && (ssec[i0]==1));
        if((strand) || (helix)) {
            segstart.push_back(i0);
            segmid.push_back(i0+((i1-i0)/2));
            segend.push_back(i1);
            seglen.push_back(i1-i0);
        }
    }

    nseg = segstart.size();
}

void Protein2::wrpdb(string fn){
    FILE *fp = openfile(fn, "w");
    wrpdb(fp);
    fclose(fp);
}

void Protein2::wrpdb(FILE *fp){
    char fmt1[] = "ATOM%7d  %-4s%3s %c%4d    %8.3f%8.3f%8.3f\n";
    for(int i=0; i<nres; i++){
        fprintf(fp, fmt1, i+1, "CA", rnam3_std[resid[i]], 'A',
                i, x[i][0], x[i][1], x[i][2]);
    }
    fprintf(fp, "TER\nEND\n");
}

double Protein2::distance2(int ia, int ib){
    return x[ia].distance2(x[ib]);
}

//void Protein2::wrbin1(std::ofstream &ofs){
int64_t Protein2::wrbin1(std::ofstream &ofs){
    //std::ofstream ofs(fn, std::ios::binary);
    const size_t n = x.size();
    const size_t n_nseg = nseg;

    ofs.write(reinterpret_cast<const char*>(&n), sizeof(n));
    ofs.write(reinterpret_cast<const char*>(&n_nseg), sizeof(n_nseg));

    //x
    for (size_t i=0; i<n; i++) {
        ofs.write(reinterpret_cast<const char*>(x[i].getx()), 3 * sizeof(double));
    }

    //x_ideal
    for (size_t i=0; i<n; i++) { //current implementation has redundant coords
        ofs.write(reinterpret_cast<const char*>(x_ideal[i].getx()), 3 * sizeof(double));
    }

    //ssec0
    ofs.write(reinterpret_cast<const char*>(ssec.data()), n*sizeof(int));
    //resid
    ofs.write(reinterpret_cast<const char*>(resid.data()), n*sizeof(int));
    // start
    ofs.write(reinterpret_cast<const char*>(segstart.data()), n_nseg*sizeof(int));
    // mid
    ofs.write(reinterpret_cast<const char*>(segmid.data()), n_nseg*sizeof(int));
    // end
    ofs.write(reinterpret_cast<const char*>(segend.data()), n_nseg*sizeof(int));
    // len
    ofs.write(reinterpret_cast<const char*>(seglen.data()), n_nseg*sizeof(int));

    return sizeof(n)+sizeof(n_nseg)+(3*2*n*sizeof(double))+(2*n*sizeof(int))+(4*n_nseg*sizeof(int));
}

void Protein2::wrbin(string fn){
    std::ofstream ofs(fn, std::ios::binary);
    wrbin1(ofs);
}

void Protein2::rdbin1(std::ifstream &ifs){
    size_t n;
    size_t n_nseg;

    size_t temp1[2];
    ifs.read(reinterpret_cast<char*>(temp1), 2*sizeof(size_t));

    memcpy(&n, temp1, sizeof(size_t));
    memcpy(&n_nseg, &temp1[1], sizeof(size_t));

	x.resize(n);
	x_ideal.resize(n);
	ssec.resize(n);
	resid.resize(n);

    segstart.resize(n_nseg);
    segmid.resize(n_nseg);
    segend.resize(n_nseg);
    seglen.resize(n_nseg);

    double temp[2*3*n];
    ifs.read(reinterpret_cast<char*>(temp), 2*3*n*sizeof(double));
    for (size_t i = 0; i < n; ++i) {
        x[i].setx(&temp[i*3]);
    }

    for (size_t i = n; i < (n+n); ++i) {
        x_ideal[i-n].setx(&temp[i*3]);
    }

    size_t total_ints = (2*n)+(4*n_nseg);
    int temp2[total_ints];

    ifs.read(reinterpret_cast<char*>(temp2), total_ints*sizeof(int));
    memcpy(ssec.data(), temp2, n*sizeof(int));
    memcpy(resid.data(), &temp2[n], n*sizeof(int));
    memcpy(segstart.data(), &temp2[2*n], n_nseg*sizeof(int));
    memcpy(segmid.data(), &temp2[(2*n)+(n_nseg)], n_nseg*sizeof(int));
    memcpy(segend.data(), &temp2[(2*n)+(2*n_nseg)], n_nseg*sizeof(int));
    memcpy(seglen.data(), &temp2[(2*n)+(3*n_nseg)], n_nseg*sizeof(int));

    nres = n;
    nseg = n_nseg;
}

void Protein2::rdbin(string fn){
    // can optimize IO by reading all together and then splitting
    std::ifstream ifs(fn, std::ios::binary);
    if(ifs.fail()) return;
    name = basename(fn.c_str());
    rdbin1(ifs);
}
