// -----------------------------------------------------------------------
// progressbar.hh: Declaration for the ProgressBar class.
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

#ifndef PROGRESSBAR_HH_
#define PROGRESSBAR_HH_

#include <iostream>
#include <string>
#include <sstream>
#include <ctime>
#include <sys/ioctl.h>


class ProgressBar           /// A class that prints out a progress 
{                           ///  bar in the form: |###           |      
public:
    ProgressBar(bool Time=false, bool Verbose=true, bool Showbar=true,
                int nlength=20, std::string ss="#");
    virtual ~ProgressBar() {}                       ///< Destructor.
    enum POS {BEG=0,END};                           ///< So that we can record where we are.

    void setShowbar (bool b) {showbar=b;}
    void setVerbose (bool b) {verbose=b;}

    void init(std::string someString, int size);    ///< Prints title, an empty bar, defines increment.
    void update(int num);                           ///< Prints correct number of hashes
    void rewind();                                  ///< Prints backspaces over bar.
    void remove();                                  ///< Overwrites bar with blanks
    void fillSpace(std::string someString);         ///< Overwrites bar with a string.
    
private:
    POS     loc;                                    ///< Are we at the start or end?
    float   stepSize;                               ///< What is the interval between hashes?
    int     length;                                 ///< What's the maximum number of hashes?
    int     numVisible;                             ///< How many hashes are there visible?
    int     stepMade;                               ///< How many steps have been completed?
    int     twidth;                                 ///< Left time string width.
    int     backs;                                  ///< Number of backspaces
    bool    ltime;                                  ///< Estimating time to finish;
    bool    showbar;                                ///< Whether to show the progress bar
    bool    verbose;                                ///< Whether to show the title and other messages.
    std::string s;                                  ///< String that fill the bar.
    struct winsize w;                               ///< Getting the size of the shell.
    clock_t start;                                  ///< Time of beginning.
    
    void defaults ();
    std::string getTimeLeft (clock_t stop);         ///< Estimating the time of the loop.
    void printBackSpaces(int n, std::ostream &str=std::cout) {for(int i=0;i<n;i++) str<<'\b';}
    void printSpaces    (int n, std::ostream &str=std::cout) {for(int i=0;i<n;i++) str<<' ';}
    void printString    (int n, std::ostream &str=std::cout) {for(int i=0;i<n;i++) str<< s;}

};


#endif
