// -----------------------------------------------------------------------
// progressbar.hh: Definitions & functions for the ProgressBar class.
// -----------------------------------------------------------------------

/*-----------------------------------------------------------------------
 This program is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation; either version 2 of the License, or (at your
 option) any later version.

 Bbarp;p is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 for more details.

 You should have received a copy of the GNU General Public License
 along with Bbarolo; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

 Correspondence concerning Bbarolo may be directed to:
    Internet email: enrico.diteodoro@unibo.it
-----------------------------------------------------------------------*/

#include <iostream>
#include <string>
#include <iomanip>
#include <ctime>
#include <unistd.h>
#include <sys/ioctl.h>
#include "progressbar.hh"
#include "utils.hh"

void ProgressBar::defaults() {
	
	length=20; 
	loc=BEG; 
	numVisible = stepMade = 0;
	s = "#";
	title = "";
	ltime = false;
    showbar = true;
	
	ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);

}

ProgressBar::ProgressBar() {
	
	defaults();
	
}


ProgressBar::ProgressBar(std::string Title, bool Time, int nlength, std::string ss) {
	
	/// This constructor enables the user to define how many string should
	/// appear. The number visible is set to 0 and the location to be at the beginning.  
	/// 
	/// \param nlength 	The new number of char to appear in the bar.
	/// \param ss 		The character to appear in the bar.
	
	defaults();
	
	length = nlength;
	s = ss.at(0);
	title = Title;
	ltime = Time;
}


void ProgressBar::init(int size) { 
	
	/// This initialises the bar to deal with a loop of a certain size.
	/// This size will imply a certain step size, dependent on the number
	/// of hashes that will be written.  A blank bar is written out as
	/// well, and we remain at the end.  
	/// 
	/// \param size 	The maximum number of iterations to be covered by the
	/// 				progress bar.

	stepSize = float(size) / float(length);
	std::cout << title << std::flush;
    if (showbar) {std::cout << "|";
        printSpaces(length);
        std::cout << "|" << std::flush;
        loc = END;
    }
		
	if (title.size()==0) twidth = 20;
	else twidth = w.ws_col-length-title.size()-8-2;
	

}


void ProgressBar::update(int num) {
	
	/// This makes sure the correct number of hashes are drawn.
	/// 
	/// Based on the number provided, as well as the stepsize, we compare
	/// the number of hashes we expect to see with the number that are
	/// there, and if they differ, the correct number are drawn. Again,
	/// we remain at the end.  
	/// 
	/// \param num 	The loop counter to be translated into the
	/// 				progress bar.
	
    if (!showbar) return;

	int numNeeded = 0;
	for(int i=0;i<length;i++)
		if(num>(i*stepSize)) numNeeded++;
	
	if(numNeeded != numVisible){
		numVisible = numNeeded;
		if(loc==END) printBackSpaces(length+2); 
		std::cout << "|"; 
		printString(numNeeded);
		printSpaces(length-numNeeded);
		std::cout << "|" << std::flush;
		loc=END;
	}
	
	if (ltime) {
	
		static clock_t first;
		clock_t second;	
		std::string timestring;
		bool flagUpdate=false;
		
		switch (stepMade) {
			case 0:
				start = clock();
				first = start;
				timestring = getTimeLeft(first);
				break;
				
			case 1:
				second = clock();
				flagUpdate=true;
				break;
			
			default:
				second = clock();
				double timestep = (double(second-first)/CLOCKS_PER_SEC);
				flagUpdate = timestep>0.99;
				break;
		}
		
		if (flagUpdate) {	
			first = second;
			timestring = getTimeLeft(second);
			std::cout << std::fixed << std::setprecision(1);	
			std::cout << std::setw(6) << std::right << float(num/(stepSize*length)*100) 
		  	  	  	  << std::setw(2) << std::right << " %"
		  	  	  	  << std::setw(twidth) << std::right << timestring << std::flush;
			printBackSpaces(twidth+8);	
		}
	}
	
	stepMade++;
}


void ProgressBar::rewind() {
	
	/// If we are at the end, we print out enough backspaces to wipe out
	/// the entire bar.  If we are not, the erasing does not need to be
	/// done.
	
	if(loc==END) printBackSpaces(length+2); 
	loc=BEG;
	std::cout << std::flush;
}


void ProgressBar::remove() {
	
	/// We first rewind() to the beginning, overwrite the bar with blank spaces, 
	/// and then rewind(). We end up at the beginning.
    if (showbar) {
        rewind();
        printSpaces(length+twidth+10);
        printBackSpaces(twidth+8);
        loc=END;
        rewind();
        std::cout << std::flush;
    }
}


void ProgressBar::fillSpace(std::string someString) {
	
	/// First remove() the bar and then write out the requested string.
	/// \param someString The string to be written over the bar area.
	
    remove();
	std::cout << someString;
	loc=END;
}


std::string ProgressBar::getTimeLeft (clock_t stop) {
	
	/// Define the time needed to end the loop. Return the time as 
	/// format hh mm ss
	
	std::string tstr=" ";
	
	double timestep = (double(stop-start)/CLOCKS_PER_SEC)/double(stepMade);
	double lefttime = timestep*(stepSize*length-stepMade);

	if (lefttime/60.>1) {
		if (lefttime/3600.>1) {
			int hours = int(lefttime/3600);
			int min = (int(lefttime)%3600)/60;
			int sec = (int(lefttime)%3600)%60;
			tstr = to_string(hours)+"h"; 
			tstr = min<10 ? tstr+"0"+to_string(min)+"m" : tstr+to_string(min)+"m";
			tstr = sec<10 ? tstr+"0"+to_string(sec)+"s" : tstr+to_string(sec)+"s";
		}
		else {
			int min = int(lefttime/60);
			int sec = int(lefttime)%60;
			tstr = to_string(min)+"m";
			tstr = sec<10 ? tstr+"0"+to_string(sec)+"s" : tstr+to_string(sec)+"s";
		}
	}
	else tstr = to_string(int(lefttime))+"s";
	
	return tstr;
}
