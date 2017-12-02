// -------------------------------------------------------------------------
// stats.hh: Definition of the Stats class, that holds statistical 
//           parameters such as mean, median, rms, madfm.
// -------------------------------------------------------------------------

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

#ifndef STATS_HH_
#define STATS_HH_

#include <iostream>
#include <cmath>


namespace Statistics
{

    /// Divide by the following correction factor to convert from MADFM to sigma (rms) estimator. 
    const float correctionFactor = 0.6744888;

    /// Multiply by the following correction factor to convert from trimmedSigma to sigma estimator. 
    const double trimToNormal = 1.17036753077;

    /// A non-templated function to do the rms-to-MADFM conversion. 
    float sigmaToMADFM(float sigma);
    
    /// A templated function to do the MADFM-to-rms conversion. 
    template <class T> float madfmToSigma(T madfm);
  
  

    ///  Template class to hold statistics for a given set of values.
    template <class Type> 
    class Stats                         
    {
    public:
        Stats() {useRobust=true; defined=false; useFDR=false;};                 /// Default constructor.
        virtual ~Stats() {};                                                    /// Default destructor.         
        Stats(const Stats<Type>& s);                                            /// Copy constructor.
        Stats<Type>& operator= (const Stats<Type>& s);                          /// Assignement operator.
        
        
        /// Obvious inline functions to access a private member of class.
        
        Type  getMean(){return mean;};
        void  setMean(Type f){mean=f;};
        Type  getStddev(){return stddev;};
        void  setStddev(Type f){stddev=f;};
        Type  getMedian(){return median;};
        void  setMedian(Type f){median=f;};
        Type  getMadfm(){return madfm;};
        void  setMadfm(Type f){madfm=f;};
        Type  getMin(){return min_val;};
        void  setMin(Type f){min_val=f;};
        Type  getMax(){return max_val;};
        void  setMax(Type f){max_val=f;};
        Type  getThreshold(){return threshold;};
        void  setThreshold(Type f){threshold=f;};
        Type  getPThreshold(){return pThreshold;};
        void  setPThreshold(Type f){pThreshold=f;};
        bool  getRobust(){return useRobust;};
        void  setRobust(bool b){useRobust=b;};
        bool  getUseFDR(){return useFDR;};
        void  setUseFDR(bool b){useFDR=b;};

        /// Specific functions.
        
        Type  getMiddle();                                      /// Return the estimator of the middle value of the data. 
        void  setMiddle(Type middle);                           /// Set the middle value of the data.               
       
        Type  getSpread();                                      /// Return the estimator of the amount of spread of the data.
        void  setSpread(Type spread);                           /// Set the amount of spread of the data.
        
        Type  getThresholdSNR();                                /// Return the threshold as a signal-to-noise ratio.
        void  setThresholdSNR(Type snr);                        /// Set the threshold in units of a signal-to-noise ratio. 
        Type  valueToSNR(Type value);                           /// Convert a value to a signal-to-noise ratio. 
        Type  snrToValue(Type snr);                         /// Convert a signal-to-noise ratio to a flux value.
        void  scaleNoise(Type scale);                           /// Scale the noise by a given factor. 
        Type  getPValue(Type value);                            /// Return the Gaussian probability of a value given the stats.  
        bool  isDetection(Type value);                          /// Is a value above the threshold? 
    
        void calculate(Type *array, long size);                 /// Calculate statistics for all elements of a data array.
        void calculate(Type *array, long size, bool *mask);     /// Calculate statistics for a subset of a data array. 
    
        template <class T> 
        friend std::ostream& operator<<(std::ostream& theStream, Stats<T> &s);  /// Overloaded operator <<.
    
    private:
        bool   defined;         ///< Are statistics defined?

        Type   mean;            ///< The mean of the values.
        Type   stddev;          ///< The standard deviation (rms).
        Type   median;          ///< The median of the values.
        Type   madfm;           ///< The median absolute deviation from the median.
        Type   max_val;         ///< Maximum value.
        Type   min_val;         ///< Minimum value.

        Type   threshold;       ///< Threshold for sigma-clipping.
        Type   pThreshold;      ///< Threshold for the FDR case (upper limit of P values that detected pixels can have).
        bool   useRobust;       ///< Use robust statistics?
        bool   useFDR;          ///< Use FDR method?

    };

}

//#include "stats.cpp"

#endif
