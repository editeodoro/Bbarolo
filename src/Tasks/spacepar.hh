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

#ifndef SPACEPAR_HH_
#define SPACEPAR_HH_

#include <iostream>
#include <cstring>
#include <Arrays/cube.hh>
#include <Tasks/galfit.hh>

using namespace Model;

template <class T>
class Spacepar : public Galfit<T>
{   
public:
    Spacepar(Cube<T> *c);
    ~Spacepar(){if (parspAll) delete parspace;}
    
    void calculate();
    void plotAll_Python();
    
private:    

    Cube<T> *parspace;  // A cube object with all parameter spaces
    bool    parspAll;   // Wheter parspace has been allocated.
    std::string p1;     // Type of parameter 1
    std::string p2;     // Type of parameter 2
    T   minp1;          // Minimum of parameter 1
    T   minp2;          // Minimum of parameter 2
    T   maxp1;          // Maximum of parameter 1
    T   maxp2;          // Maximum of parameter 2
    T   delp1;          // Step size for parameter 1
    T   delp2;          // Step size for parameter 2
    
    void plotmin_Gnuplot(int n);
    void setOutSpacepar(int p1_nsteps, int p2_nsteps);
};




#endif
