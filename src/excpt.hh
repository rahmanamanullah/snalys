#ifndef INCLUDE_EXCPT_HH
#define INCLUDE_EXCPT_HH 1

#include <string>

using namespace std;

class Excpt {
public:
  string errclass,errmethod,errmsg;
  
  Excpt(const string &c,const string &m,const string &s) {
    errclass=c;
    errmethod=m;
    errmsg=s;
  }
};

#endif
