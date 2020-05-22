/////////////////////////////////////////////
// bbarolo.cpp: Main source file of BBarolo
/////////////////////////////////////////////

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
#include <sys/time.h>
#include <iomanip>
#include <bbarolo.hh>
#include <Arrays/param.hh>


/*
#include<signal.h>
struct sigaction osa;
void action_sigint(int sig_no)
{
    printf("\nI tap SIGINT and returns back \n");
    sigaction(SIGINT,&osa,NULL);
    kill(0,SIGINT);
}
*/
int main (int argc, char *argv[]) {

//    struct sigaction act;
//    act.sa_handler = action_sigint;
//    sigemptyset(&act.sa_mask);
//    act.sa_flags = 0;
//    sigaction(SIGINT, &act, 0);
    
    struct timeval begin, end;
    gettimeofday(&begin, NULL);

    Param *par = new Param;
    if (!par->getopts(argc, argv)) return EXIT_FAILURE;
    bool verbose = par->isVerbose();
    if (verbose) {
        welcomeMessage();
        std::cout << *par;
    }

    for (int im=0; im<par->getListSize(); im++) {

        if (par->getListSize()>1 && verbose) {
            std::cout << setfill('_') << std::endl;
            std::cout << setw(70) << "" << std::endl << std::endl;
            std::string s = "Working on "+ par->getImage(im)+" ";
            std::cout << setfill(' ') << right << setw(70) << s;
            std::cout << std::endl << left;
            std::cout << std::endl << " File "<< im+1 << " of " 
                      << par->getListSize()<< std::endl << std::endl;
        }

        par->setImageFile(par->getImage(im));

        if (!BBcore(par)) {
            if(par->getListSize()-im>1) std::cout << "Skipping to next file...\n";
            else {std::cout << "Exiting ...\n\n"; return EXIT_FAILURE;}
        }

    }

    delete par;
    
    gettimeofday(&end, NULL);
    double time = (end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0);
    if (verbose)
        std::cout << "\nExecution time: " << int(time/60) << " min and " << int(time)%60 << " sec.\n";

    return EXIT_SUCCESS;

}
