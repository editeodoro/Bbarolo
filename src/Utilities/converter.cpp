// ----------------------------------------------------------------
//	converter.cpp: Definitions of functions for the UNITS class
// ----------------------------------------------------------------

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
#include <cmath>
#include "converter.hh"

double UNITS::ConvFactor() {
	
	bool InIsLinear=false, OutIsLinear=false;
	for (int i=ANG;i<=MPC;i++) {
		InIsLinear  = InIsLinear  || in==i;
		OutIsLinear = OutIsLinear || out==i;
	}
	if (InIsLinear && OutIsLinear) return LinearFactor();
	
	bool InIsTime=false, OutIsTime=false;
	for (int i=S;i<=GYR;i++) {
		InIsTime  = InIsTime  || in==i;
		OutIsTime = OutIsTime || out==i;
	}
	if (InIsTime && OutIsTime) return TimeFactor();
	
	bool InIsMass=false, OutIsMass=false;
	for (int i=GR;i<=MSUN;i++) {
		InIsMass  = InIsMass  || in==i;
		OutIsMass = OutIsMass || out==i;
	}
	if (InIsMass && OutIsMass) return MassFactor();
	
	bool InIsFreq=false, OutIsFreq=false;
	for (int i=HZ;i<=GHZ;i++) {
		InIsFreq  = InIsFreq  || in==i;
		OutIsFreq = OutIsFreq || out==i;
	}
	if (InIsFreq && OutIsFreq) return FrequencyFactor();
	
	bool InIsAngle=false, OutIsAngle=false;
	for (int i=DEG;i<=RAD;i++) {
		InIsAngle  = InIsAngle  || in==i;
		OutIsAngle = OutIsAngle || out==i;
	}
	if (InIsAngle && OutIsAngle) return AngleFactor();
	
	bool InIsVel=false, OutIsVel=false;
	for (int i=CM_S;i<=KM_H;i++) {
		InIsVel  = InIsVel  || in==i;
		OutIsVel = OutIsVel || out==i;
	}
	if (InIsVel && OutIsVel) {
		Unit NumIn=M, NumOut=M, DenIn=S, DenOut=S;
		if (in==CM_S) NumIn = CM; 
		else if (in==M_S) NumIn = M;
		else if (in==KM_S || in==KM_H) NumIn = KM;	
		if (out==CM_S) NumOut = CM; 
		else if (out==M_S) NumOut = M;
		else if (out==KM_S || out==KM_H) NumOut = KM;
		if (in==CM_S || in==M_S || in==KM_S) DenIn = S;
		else if (in==KM_H) DenIn = H;
		if (out==CM_S || out==M_S || out==KM_S) DenOut = S;
		else if (out==KM_H) DenOut = H;
		UNITS NUM(NumIn,NumOut), DENOM(DenIn,DenOut);
		return NUM.factor/DENOM.factor;
	}

	std::cout << "Unknown conversion!!" << std::endl;
	
	return -1;
}


double UNITS::LinearFactor () {
	
	double factor = 1;
	
	/// First convert all to meters
	switch (in) {
		case ANG:
			factor = 1E-10;
			break;
		case NM:
			factor = 1E-09;
			break;
		case MUM:
			factor = 1E-06;
			break;
		case MM:
			factor = 1E-03;
			break;
		case CM:
			factor = 1E-02;
			break;
		case M:
			factor = 1;
			break;
		case KM:
			factor = 1E+03;
			break;
		case UA:
			factor = 1.49597870700E11;
			break;
		case PC:
			factor = 3.08567758E16;
			break;
		case KPC:
			factor = 3.08567758E19;
			break;
		case MPC:
			factor = 3.08567758E22;
			break;				
		default:
			ErrorMessage();
			return factor;
	}
	
	
	switch (out) {
		case ANG:
			factor *= 1E+10;
			break;
		case NM:
			factor *= 1E+09;
			break;
		case MUM:
			factor *= 1E+06;
			break;
		case MM:
			factor *= 1E+03;
			break;
		case CM:
			factor *= 1E+02;
			break;
		case M:
			factor *= 1;
			break;
		case KM:
			factor *= 1E-03;
			break;
		case UA:
			factor *= 6.68458712E-12;
			break;
		case PC:
			factor *= 3.24077929E-17;
			break;
		case KPC:
			factor *= 3.24077929E-20;
			break;
		case MPC:
			factor *= 3.24077929E-23;
			break;				
		default:
			ErrorMessage();
			factor=1;
			return factor;
	}
	
	return factor;
}


double UNITS::TimeFactor() {
	 
	double factor = 1;
	 
	// First convert all to seconds.
	switch (in) {
		case S:
			factor = 1;
			break;
		case H:
			factor = 3600;
			break;
		case D:
			factor = 86400;
			break;
		case YR:
			factor = 3.1556926E+07;
			break;
		case KYR:
			factor = 3.1556926E+10;
			break;
		case MYR:
			factor = 3.1556926E+13;
			break;
		case GYR:
			factor = 3.1556926E+16;
			break;
		default:
			ErrorMessage();
			return factor;		 
	}
	 
	switch (out) {
		case S:
			factor *= 1;
			break;
		case H:
			factor *= 0.000277777778;
			break;
		case D:
			factor *= 1.15740741E-05;
			break;
		case YR:
			factor *= 3.16887646E-08;
			break;
		case KYR:
			factor *= 3.16887646E-11;
			break;
		case MYR:
			factor *= 3.16887646E-14;
			break;
		case GYR:
			factor *= 3.16887646E-17;
			break;
		default:
			ErrorMessage();
			factor = 1;
			return factor;		 
	}
	
	return factor;
}


double UNITS::MassFactor() {
	 
	double factor = 1;
	 
	// First convert all to KG.
	switch (in) {
		case GR:
			factor = 1E-03;
			break;
		case KG:
			factor = 1;
			break;
		case MSUN:
			factor = 1.9891E30;
			break;
		default:
			ErrorMessage();
			return factor;		 
	}
	 
	switch (out) {
		case GR:
			factor *= 1E+03;
			break;
		case KG:
			factor *= 1;
			break;
		case MSUN:
			factor *= 5.02739933E-31;
			break;
		default:
			ErrorMessage();
			factor = 1;
			return factor;		 
	}
	
	return factor;
}


double UNITS::FrequencyFactor() {
	 
	double factor = 1;
	 
	// First convert all to GHZ.
	switch (in) {
		case HZ:
			factor = 1E-09;
			break;
		case KHZ:
			factor = 1E-06;
			break;
		case MHZ:
			factor = 1E-03;
			break;
		case GHZ:
			factor = 1;
			break;
		default:
			ErrorMessage();
			return factor;		 
	}
	 
	switch (out) {
		case HZ:
			factor *= 1E+09;
			break;
		case KHZ:
			factor *= 1E+06;
			break;
		case MHZ:
			factor *= 1E+03;
			break;
		case GHZ:
			factor *= 1;
			break;
		default:
			ErrorMessage();
			factor = 1;
			return factor;		 
	}
	
	return factor;
}


double UNITS::AngleFactor() {
	 
	double factor = 1;
	 
	// First convert all to DEGREE.
	switch (in) {
		case DEG:
			factor = 1;
			break;
		case AMIN:
			factor = 1.66666666667E-02;
			break;
		case ASEC:
			factor = 2.77777777778E-04;
			break;
		case RAD:
			factor = 180./M_PI;
			break;
		default:
			ErrorMessage();
			return factor;		 
	}
	 
	switch (out) {
		case DEG:
			factor *= 1;
			break;
		case AMIN:
			factor *= 60;
			break;
		case ASEC:
			factor *= 3600;
			break;
		case RAD:
			factor *= M_PI/180.;;
			break;
		default:
			ErrorMessage();
			factor = 1;
			return factor;		 
	}
	
	return factor;
}


void UNITS::ErrorMessage () {
	std::cout << "Conversion error: unknown unit " << std::endl; 
}


