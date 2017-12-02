//------------------------------------------------------------------
//  converter.hh: Conversion between physical units
//------------------------------------------------------------------

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

#ifndef CONVERTER_HH_
#define CONVERTER_HH_

enum Unit {
    ANG, NM, MUM, MM, CM, M, KM,                        // Linear units  
    UA, PC, KPC, MPC,                                   // Linear astronomical units
    S, H, D, YR, KYR, MYR, GYR,                     // Time units
    GR, KG, MSUN,                                       // Mass units
    HZ, KHZ, MHZ, GHZ,                                  // Frequency units
    DEG, AMIN, ASEC, RAD,                               // Angle units
    CM_S, M_S, KM_S, KM_H,                              // Velocity units
};


class UNITS 
{
public:
    Unit in;
    Unit out;
    double factor;
    
    UNITS(Unit IN, Unit OUT) : in(IN), out(OUT) {setFactor();}
    ~UNITS() {};
    
    void setFactor() {factor=ConvFactor();}
    
private:
    double  ConvFactor();
    double LinearFactor ();
    double TimeFactor ();
    double MassFactor ();
    double FrequencyFactor ();
    double AngleFactor ();
    void    ErrorMessage();
};

// Convert the value from Units in to Units out
template <class T> double ConvertUnits (UNITS *units, T value) {
    return value*units->factor;
}

template <class T> double ConvertUnits (Unit in, Unit out, T value) {
    return ConvertUnits(new UNITS(in,out),value);
}


// Get the conversion factor between Units in and out
template <class T> double ConversionFactor (Unit in, Unit out) {
    UNITS u(in, out);
    return u.factor;
}

#endif
