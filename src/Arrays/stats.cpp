// -----------------------------------------------------------------------
// stats.cpp: Member functions for the templated Stats class.
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
    Internet email: enrico.diteodoro@unibo.it
-----------------------------------------------------------------------*/

#include <iostream>
#include <Arrays/stats.hh>
#include <Utilities/utils.hh>


namespace Statistics
{

	float sigmaToMADFM(float sigma){
		return float(sigma)*correctionFactor;
	}


	template <class T> 
	float madfmToSigma(T madfm){
		return float(madfm)/correctionFactor;
	}
	template float madfmToSigma<short>(short);
	template float madfmToSigma<int>(int);
	template float madfmToSigma<long>(long);
	template float madfmToSigma<float>(float);
	template float madfmToSigma<double>(double);


//-----------------------------------------------------------------------------


	template <class Type> 
	Stats<Type>::Stats(const Stats<Type>& s) {
		
		operator=(s);
	}
	template Stats<short>::Stats(const Stats<short>&);
	template Stats<int>::Stats(const Stats<int>&);
	template Stats<long>::Stats(const Stats<long>&);
	template Stats<float>::Stats(const Stats<float>&);
	template Stats<double>::Stats(const Stats<double>&);



	template <class Type> 
	Stats<Type>& Stats<Type>::operator= (const Stats<Type>& s) {

		if(this == &s) return *this;
		this->defined    = s.defined;
		this->mean       = s.mean;
		this->stddev     = s.stddev;
		this->median     = s.median;
		this->madfm      = s.madfm;
		this->max_val    = s.max_val;
		this->min_val    = s.min_val;
		this->threshold  = s.threshold;
		this->pThreshold = s.pThreshold;
		this->useRobust  = s.useRobust;
		this->useFDR     = s.useFDR;

		return *this;
	}
	template Stats<short>& Stats<short>::operator= (const Stats<short>&);
	template Stats<int>& Stats<int>::operator= (const Stats<int>&);
	template Stats<long>& Stats<long>::operator= (const Stats<long>&);
	template Stats<float>& Stats<float>::operator= (const Stats<float>&);
	template Stats<double>& Stats<double>::operator= (const Stats<double>&);



	template <class Type> 
	Type Stats<Type>::getThresholdSNR() {
			return (threshold - getMiddle())/getSpread();	/// The SNR is defined in terms of excess over the middle 
	}														/// estimator in units of the spread estimator.	
	template short Stats<short>::getThresholdSNR();
	template int Stats<int>::getThresholdSNR();
	template long Stats<long>::getThresholdSNR();
	template float Stats<float>::getThresholdSNR();
	template double Stats<double>::getThresholdSNR();



	template <class Type> 
	void  Stats<Type>::setThresholdSNR(Type snr) {
		threshold = getMiddle() + snr*getSpread();
	}
	template void Stats<short>::setThresholdSNR(short);
	template void Stats<int>::setThresholdSNR(int);
	template void Stats<long>::setThresholdSNR(long);
	template void Stats<float>::setThresholdSNR(float);
	template void Stats<double>::setThresholdSNR(double);


  
	template <class Type> 
	Type Stats<Type>::valueToSNR(Type value) {
		return (value - getMiddle())/getSpread();
	}
	template short Stats<short>::valueToSNR(short);
	template int Stats<int>::valueToSNR(int);
	template long Stats<long>::valueToSNR(long);
	template float Stats<float>::valueToSNR(float);
	template double Stats<double>::valueToSNR(double);
 
 
  
	template <class Type> 
	Type Stats<Type>::snrToValue(Type snr) {
		return (snr * getSpread() + getMiddle());
	}
	template short Stats<short>::snrToValue(short);
	template int Stats<int>::snrToValue(int);
	template long Stats<long>::snrToValue(long);
	template float Stats<float>::snrToValue(float);
	template double Stats<double>::snrToValue(double);
 
 
    
	template <class Type> 
	void Stats<Type>::setMiddle(Type middle) {
   
  /// The middle value is determined by the Stats::useRobust flag 
  /// It will be either the median (if true), or the mean (if false).
  
    if(useRobust) median = Type(middle); 
    else mean = middle;
  }
	template void Stats<short>::setMiddle(short);
	template void Stats<int>::setMiddle(int);
	template void Stats<long>::setMiddle(long);
	template void Stats<float>::setMiddle(float);
	template void Stats<double>::setMiddle(double);


    
	template <class Type> 
	Type Stats<Type>::getMiddle() {
		if(useRobust) return median; 
		else return mean;
	}
	template short Stats<short>::getMiddle();
	template int Stats<int>::getMiddle();
	template long Stats<long>::getMiddle();
	template float Stats<float>::getMiddle();
	template double Stats<double>::getMiddle();
 
    
	template <class Type> 
	void Stats<Type>::setSpread(Type spread) {
  
  /// The spread value is set according to the Stats::useRobust flag --
  /// It will be either the madfm (if true), or the rms (if false). 
  /// If robust, the spread value will be converted to a madfm from an 
  /// equivalent rms under the assumption of Gaussianity, using the 
  /// Statistics::sigmaToMADFM function.
  
		if(useRobust) madfm = Type(sigmaToMADFM(spread)); 
		else stddev = spread;
  }
	template void Stats<short>::setSpread(short);
	template void Stats<int>::setSpread(int);
	template void Stats<long>::setSpread(long);
	template void Stats<float>::setSpread(float);
	template void Stats<double>::setSpread(double);
  
  
     
	template <class Type> 
	Type Stats<Type>::getSpread() {

  /// The spread value returned is determined by the Stats::useRobust flag. 
  /// It will be either the madfm (if true), or the rms (if false). 
  /// If robust, the madfm will be converted to an equivalent rms under the 
  /// assumption of Gaussianity, using the Statistics::madfmToSigma function.
		
		if(useRobust) return madfmToSigma(madfm); 
		else return stddev;
	}
	template short Stats<short>::getSpread();
	template int Stats<int>::getSpread();
	template long Stats<long>::getSpread();
	template float Stats<float>::getSpread();
	template double Stats<double>::getSpread();
  
  
    
	template <class Type> 
	void  Stats<Type>::scaleNoise(Type scale) {
	  
  ///  Multiply the noise parameters (stddev & madfm) by a given
  ///  factor, and adjust the threshold.
   
		Type snr = (threshold - getMiddle())/getSpread();    
		madfm  = Type(madfm*scale);
		stddev *= scale;
		threshold = getMiddle() + snr*getSpread();
	}
	template void Stats<short>::scaleNoise(short);
	template void Stats<int>::scaleNoise(int);
	template void Stats<long>::scaleNoise(long);
	template void Stats<float>::scaleNoise(float);
	template void Stats<double>::scaleNoise(double);



	template <class Type> 
	Type Stats<Type>::getPValue(Type value) {

    /// Get the probability, under the assumption of normality, of a
    /// value occuring.  
    
		Type zStat = (value-getMiddle())/getSpread();
		return 0.5*erfc(zStat/M_SQRT2);
	}
	template short Stats<short>::getPValue(short);
	template int Stats<int>::getPValue(int);
	template long Stats<long>::getPValue(long);
	template float Stats<float>::getPValue(float);
	template double Stats<double>::getPValue(double);


    
	template <class Type> 
	bool Stats<Type>::isDetection(Type value) {

  /// Compares the value given to the correct threshold, depending on
  /// the value of the Stats::useFDR flag. 
    
		if(useFDR) return (getPValue(value) < pThreshold);
		else return (value > threshold);
	}
	template bool Stats<short>::isDetection(short);
	template bool Stats<int>::isDetection(int);
	template bool Stats<long>::isDetection(long);
	template bool Stats<float>::isDetection(float);
	template bool Stats<double>::isDetection(double);



	template <class Type> 
	void Stats<Type>::calculate(Type *array, long size) {

    /// Calculate all four statistics for all elements of a given array.
    /// 
    /// \param array 	The input data array.
    /// \param size 	The length of the input array
    
		findAllStats(array,size,mean,stddev,median,madfm,max_val,min_val);   	
		defined = true;
	}
	template void Stats<short>::calculate(short *, long);
	template void Stats<int>::calculate(int *, long);
	template void Stats<long>::calculate(long *, long);
	template void Stats<float>::calculate(float *, long);
	template void Stats<double>::calculate(double *, long);



	template <class Type> 
	void Stats<Type>::calculate(Type *array, long size, bool *mask) {

  /// Calculate all four statistics for a subset of a given array. 
  /// The subset is defined by an array of bool variables.  
  /// 
  /// \param array 		The input data array.
  /// \param size	 	The length of the input array
  /// \param mask 		An array of the same length that says whether to
  /// 					include each member of the array in the calculations.

        findAllStats(array, size, mask, mean, stddev, median, madfm,max_val,min_val);
		defined = true;
	}
	template void Stats<short>::calculate(short *, long, bool *);
	template void Stats<int>::calculate(int *, long, bool *);
	template void Stats<long>::calculate(long *, long, bool *);
	template void Stats<float>::calculate(float *, long, bool *);
	template void Stats<double>::calculate(double *, long, bool *);



	template <class Type> 
	std::ostream& operator<< (std::ostream& theStream, Stats<Type> &s) {
		
    /// Prints out the four key statistics to the requested stream.
    
    theStream << "Mean   = "   << s.mean   << "\t"
			  << "Std.Dev. = " << s.stddev << "\n"
			  << "Median = "   << s.median << "\t"
			  << "MADFM    = " << s.madfm  << " (= " << madfmToSigma(s.madfm) << " as std.dev.)\n";
    return theStream;
	}
	template std::ostream& operator<<<short> (std::ostream& theStream, Stats<short> &s);
	template std::ostream& operator<<<int> (std::ostream& theStream, Stats<int> &s);
	template std::ostream& operator<<<long> (std::ostream& theStream, Stats<long> &s);
	template std::ostream& operator<<<float> (std::ostream& theStream, Stats<float> &s);
	template std::ostream& operator<<<double> (std::ostream& theStream, Stats<double> &s);


}
