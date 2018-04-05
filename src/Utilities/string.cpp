// -----------------------------------------------------------------------
// string.cpp: General utility functions for manipulating strings 
//                    and input flags/parameters
// -----------------------------------------------------------------------

/*-----------------------------------------------------------------------
 This program is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation; either version 2 of the License, or (at your
 option) any later version.

 BBarolo is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 for more details.

 You should have received a copy of the GNU General Public License
 along with BBarolo; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

 Correspondence concerning BBarolo may be directed to:
    Internet email: enrico.diteodoro@gmail.com
-----------------------------------------------------------------------*/

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <ctime>
#include <random>
#include <Utilities/utils.hh>


std::string makelower(std::string s) {

    std::string out = "";
    for( size_t i=0; i<s.size(); ++i ) {
        out += tolower(s[i]);
    }
    return out;
}


std::string makeupper(std::string s) {

    std::string out = "";
    for( size_t i=0; i<s.size(); ++i ) {
        out += toupper(s[i]);
    }
    return out;
}

std::string stringize(bool b) {
    
  /// Convert a bool variable to the textual equivalent. 
  /// \return A std::string with the english equivalent of the bool.

    std::string output;
    if(b) output="true";
    else output="false";
    return output;
}


bool boolify(std::string s) {
  
  ///  Convert a std::string to a bool variable: 
  ///  "1" and "true" get converted to true;
  ///  "0" and "false" (and anything else) get converted to false.
  ///  \return The bool equivalent of the string.

    if((s=="1") || (makelower(s)=="true")) return true;
    else if((s=="0") || (makelower(s)=="false")) return false;
    else return false;
}


std::string readSval(std::stringstream& ss) {
 
  std::string val; 
  ss >> val; 
  return val;
}


std::string readFilename(std::stringstream& ss) {
  
    std::string val;
    getline(ss,val);
    return deblank(val);
}


bool readFlag(std::stringstream& ss) {
    
    std::string val; 
    ss >> val; 
    return boolify(val);
}


float readFval(std::stringstream& ss) {
  
    float val; 
    ss >> val; 
    return val;
}


double readDval(std::stringstream& ss) {
  
    double val; 
    ss >> val; 
    return val;
}


int readIval(std::stringstream& ss) {
 
  int val; 
  ss >> val; 
  return val;
}


template <class T>
void readArray(std::stringstream& ss, T *val, int n) {
 
  for (int i=0; i<n; i++) ss >> val[i]; 
}
template void readArray(std::stringstream&, int*, int);
template void readArray(std::stringstream&, long*, int);
template void readArray(std::stringstream&, float*, int);
template void readArray(std::stringstream&, double*, int);


template <class T>
std::vector<T> readVec(std::stringstream& ss) {
  T value; 
  std::vector<T> val;
  while (ss >> value) val.push_back(value);
  return val;
}
template std::vector<int> readVec(std::stringstream&);
template std::vector<long> readVec(std::stringstream&);
template std::vector<float> readVec(std::stringstream&);
template std::vector<double> readVec(std::stringstream&);


std::string removeLeadingBlanks(std::string s) {
  
 /// All blank spaces from the start of the string to the first
 /// non-blank-space character are deleted.
    
    int i=0;
    while(isspace(s[i])){
        i++;
    }
    std::string newstring="";
    for(unsigned int j=i;j<s.size();j++) newstring += s[j];
    return newstring;
}


std::string deblank(std::string s) {
  
  /// All blank spaces from the start of the string to the first
  /// non-blank-space character, and from the last non-blank-space
  /// character to the end are deleted.
   
    int beg=0;
    while(isspace(s[beg])){
        beg++;
    }
    int end=s.size()-1;
    while(isspace(s[end])){
        end--;
    }
    std::string newstring;
    for(int j=beg;j<=end;j++) newstring += s[j];
    return newstring;

}


std::string deblankAll (std::string s) {
        
    // Delete all blank spaces from the string s.
    
    std::string newstring;
    for (uint i=0; i<s.size();i++) 
        if (!isspace(s[i])) newstring += s[i];
    return newstring;
    
}


template <class T> std::string to_string (const T& t, int precision) {
    
  /// Convert another type T to string.
  
    std::stringstream ss;
    ss << std::fixed;
    if (precision!=-1) ss << std::setprecision(precision);
    ss << t;
    return ss.str();
}
template std::string to_string (const short&,int);
template std::string to_string (const int&,int);
template std::string to_string (const long&,int);
template std::string to_string (const float&,int);
template std::string to_string (const double&,int);



void printBackSpace(std::ostream &stream, int num) {
  
  for(int i=0;i<num;i++) stream << '\b';
}


void printBackSpace(int num) {
    
    printBackSpace(std::cout,num);
}


void printSpace(std::ostream &stream, int num) {
    
    for(int i=0;i<num;i++) stream << ' '; 
}


void printSpace(int num) {
    
    printSpace(std::cout,num);
}


void printHash(std::ostream &stream, int num) {
    
    for(int i=0;i<num;i++) stream << '#'; 
}


void printHash(int num) {
 
    printHash(std::cout,num);
}

void checkHome(std::string &s) {
    
    if (s[0]=='~') {
        std::string home = std::getenv("HOME");
        s.erase(0,1);
        s.insert(0, home);
    }    
    /*
    if (s.find("./")==0) {
        std::string path = std::getenv("PWD");
        s.erase(0,1);
        s.insert(0, path);
    }
    */
    
}


std::string randomAdjective (int type) {
    
    // Return a random good (type=1) or bad (type=2) adjective
    
    std::vector<std::string> good = 
        {"trustful","awesome","stunning","creative","flabbergasting","fabulous","gorgeous",
         "unbelievable","extraordinary","breathtaking","astonishing","stonking","brilliant",
         "marvelous","ungodly","incredible","wondrous","magnificent","glorious","splendiferous",
         "wonderful","phantasmagoric","phenomenal","excellent","exceptional","refulgent",
         "ambitious","exuberant","frank","witty","amiable","fearless","honest","ineffable",
         "arcadian","egregious","didactic","efficacious","judicious","propitious","sagacious",
         "zealous","flamboyant"}; 
    
    std::vector<std::string> bad = 
        {"flummoxed","cranky","pernicious","modest","shameful","wobbling","knackered","flippant",
         "wonky","bellicose","caustic","calamitous","crapulous","dowdy","execrable","fastidious",
         "guileless","hubristic","insidious","insolent","irksome","mendacious","meretricious",
         "noxious","obtuse","recalcitrant","risible","strident","wheedling","withering","pauciloquent",
         "preposterous","brazenly"};
    
    static auto const seed = std::random_device()();
    static std::mt19937 generator(seed);
    
    if (type==1) {
        std::uniform_int_distribution<int> distr(0,good.size()-1);
        return good[distr(generator)];
    }
    else {
        std::uniform_int_distribution<int> distr(0,bad.size()-1);
        return bad[distr(generator)];
    }
}
