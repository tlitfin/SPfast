#ifndef _SP_MISC
#define _SP_MISC
#include <algorithm>
#include <cassert>
#include <vector>
#include <set>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <map>
#include <ctime>
using namespace std;

extern int DEBUG;
#define die0() do {fprintf(stderr, "\nERROR ON LINE %d OF FILE %s!\n\n", __LINE__, __FILE__); exit(1);} while(0)
#define die1(...) do {warn(__VA_ARGS__); die0();} while(0)
#define die(...) do {warn(__VA_ARGS__); exit(0);} while(0)
template <class T>
inline void build_2D_array(vector< vector<T> > & a, int x, int y){
    a.resize(x, vector<T>(y));
}
template <class T>
inline double dot_product(int n, T r1, T r2){
    double value=0;
    for(int i=0;i<n;i++) value += r1[i]*r2[i];
    return value;
}
template <class T>
inline double dot_product(T r1, T r2){
    return dot_product(3, r1, r2);
}
inline bool endswith(const string s1, const string s2){
    if(s1.size()>=s2.size() && s1.substr(s1.size()-s2.size(), s2.size())==s2) return true;
    return false;
}
inline FILE *openfile(const string file, const string action="r"){
    FILE *fp = NULL;
    if(file == "--"){
        if(action[0] == 'r') return stdin;
        else if(action[0] == 'w') return stdout;
    }
    if(file=="STDIN" && action[0]=='r') return stdin;
    if(file=="STDOUT" && action[0]=='w') return stdout;
    if(action == "r"){
        if(endswith(file, ".gz")) {
            string cmd1 = "gunzip -c " + file;
            return popen(cmd1.c_str(), "r");
        }
    }
//
    fp = fopen(file.c_str(), action.c_str());
    //fp = fopen(file.c_str(), "r");
    if(fp == NULL) fprintf(stderr, "Warning, Fail to open file: %s\n", file.c_str());
    return fp;
}
template <class T>
inline int maxloc(int n, T dat){
    int imax = 0;
    for(int i=1; i<n; i++){
        if(dat[i] > dat[imax]) imax = i;
    }
    return imax;
}
template <class T>
inline int str2dat(const string &line, T *dat){
   istringstream iss(line);
    int n=0;
   while(! (iss>>dat[n]).fail() ) n++;
    return n;
}
template <class T>
inline bool str2dat(const string &line, const int n, T* dat){
   istringstream iss(line);
    for(int i=0; i<n; i++){
    if( (iss>>dat[i]).fail() )  return false;
    }
    return true;
}
template <class T>
inline int str2dat(const string &line, vector<T> &dat){
   istringstream iss(line);
    T dt; int n=0, m=dat.size();
   while(! (iss>>dt).fail() ){
        if(n < m) dat[n] = dt;
        else dat.push_back(dt);
        n ++;
    }
    return n;
}
inline void warn(const char *fmt, ...){
    va_list ptr; va_start(ptr, fmt);
    vfprintf(stderr, fmt, ptr);
    va_end(ptr);
    fprintf(stderr, "\n");
}
template <class T>
inline T powi(T x, int n){
    if(n == 0) return 1.;
    if(n < 0) return 1. / powi(x, -n);
    T w = x, y = 1;
    if (n & 1) y = x;
    n >>= 1;
    while(n) {
        w = w * w;
        if( n & 1 ) y *= w;
        n >>= 1;
    }
    return y;
}
inline char* chomp(char *str){
    int nt = strlen(str);
    if(str[nt-1] == '\n') str[nt-1] = '\0';
    return str;
}
inline bool file_existed(string fn){
    FILE *fp = fopen(fn.c_str(), "r");
    if(fp == NULL) return false;
    fclose(fp); return true;
}
#endif
