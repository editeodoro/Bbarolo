// -----------------------------------------------------------------------
// statistics.ccp: Functions to calculate statistical parameters
//            	   (mean, median, rms, madfm, min, max).
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
#include <algorithm>
#include <cmath>
#include <vector>
#include "utils.hh"


template <class T> 
T absval(T value) {
  
  /// Type-independent way of getting the absolute value.
	if( value > 0) return value;
	else return 0-value;
}
template short absval<short>(short);
template int absval<int>(int);
template long absval<long>(long);
template float absval<float>(float);
template double absval<double>(double);



template <class T> 
void findMinMax(const T *array, const size_t size, T &min, T &max) {

  /// A function to find the minimum and maximum values of a set of numbers.
  ///
  /// \param array 	The array of data values. Type independent.
  /// \param size 	The length of the array
  /// \param min 	The returned value of the minimum value in the array.
  /// \param max 	The returned value of the maximum value in the array.
  
	min = max = array[0];
	for(size_t i=1;i<size;i++) {
        if (array[i]!=array[i]) continue;
		if(array[i]<min) min=array[i];
		if(array[i]>max) max=array[i];
	}
}
template void findMinMax<short>(const short *, const size_t, short &, short &);
template void findMinMax<int>(const int *, const size_t, int &, int &);
template void findMinMax<long>(const long *, const size_t, long &, long &);
template void findMinMax<float>(const float *, const size_t, float &, float &);
template void findMinMax<double>(const double *, const size_t, double &, double &);



template <class T> 
void findAllStats(T *array, size_t size, T &mean, T &stddev, T &median, T &madfm) {

  /// Find the mean,rms (or standard deviation), median AND madfm of an
  /// array of numbers. Type independent.
  /// 
  /// \param array 	The array of numbers.
  /// \param size 	The length of the array.
  /// \param mean 	The mean value of the array, returned as a float.
  /// \param stddev The rms or standard deviation of the array,
  /// 				returned as a float.
  /// \param median The median value of the array, returned as the same
  /// 				type as the array.
  /// \param madfm 	The median absolute deviation from the median value
  /// 			   	of the array, returned as the same type as the array.

	T maxx, minn;
	findAllStats(array,size,mean,stddev,median,madfm,maxx,minn);

}
template void findAllStats<short>(short *, size_t, short&, short&, short&, short&);
template void findAllStats<int>(int *, size_t, int&, int&, int&, int&);
template void findAllStats<long>(long *, size_t, long&, long&, long&, long&);
template void findAllStats<float>(float *, size_t, float&, float&, float&, float&);
template void findAllStats<double>(double *, size_t, double&, double&, double&, double&);



template <class T> 
void findAllStats(T *array, size_t size, T &mean, T &stddev, T &median, T &madfm, T &maxx, T &minn) {

  /// Find the mean,rms (or standard deviation), median AND madfm of an
  /// array of numbers. Type independent.
  /// 
  /// \param array 	The array of numbers.
  /// \param size 	The length of the array.
  /// \param mean 	The mean value of the array, returned as a float.
  /// \param stddev The rms or standard deviation of the array,
  /// 				returned as a float.
  /// \param median The median value of the array, returned as the same
  /// 				type as the array.
  /// \param madfm 	The median absolute deviation from the median value
  /// 			   	of the array, returned as the same type as the array.

	if(size==0){
		std::cout << "Error in findAllStats: zero sized array!\n";
		return;
	}

	T *newarray = new T[size];

	minn = maxx = array[0];
	for(size_t i=0;i<size;i++) {
		if(array[i]<minn) minn=array[i];
		if(array[i]>maxx) maxx=array[i];
		newarray[i] = array[i];
	}

	mean = findMean(newarray,size);
	stddev = findStddev(newarray,size);
	median = findMedian(newarray,size,true);
	madfm = findMADFM(newarray,size,median,true);

	delete [] newarray;

}
template void findAllStats<short>(short *, size_t, short&, short&, short&, short&, short&, short&);
template void findAllStats<int>(int *, size_t, int&, int&, int&, int&, int&, int&);
template void findAllStats<long>(long *, size_t, long&, long&, long&, long&, long&, long&);
template void findAllStats<float>(float *, size_t, float&, float&, float&, float&, float&, float&);
template void findAllStats<double>(double *, size_t, double&, double&, double&, double&, double&, double&);



template <class T> 
void findAllStats(T *array, size_t size, bool *mask, T &mean, T &stddev, T &median, T &madfm, T &maxx, T &minn) {

  /// Find the mean,rms (or standard deviation), median AND madfm of a
  /// subset of an array of numbers. Type independent. The subset is
  /// defined by an array of bool variables. Type independent.
  /// 
  /// \param array 	The array of numbers.
  /// \param size 	The length of the array.
  /// \param mask 	An array of the same length that says whether to
  /// 			  	include each member of the array in the calculations. 
  ///			  	Only look at values where mask=true.
  /// \param mean 	The mean value of the array, returned as a float.
  /// \param stddev The rms or standard deviation of the array,
  /// 				returned as a float.
  /// \param median The median value of the array, returned as the same
  /// 				type as the array.
  /// \param madfm 	The median absolute deviation from the median value
  /// 			   	of the array, returned as the same type as the array.

	int goodSize=0;
	for(size_t i=0;i<size;i++) if(mask[i]) goodSize++;
	if(goodSize==0){
		std::cout << "Error in findAllStats: no good values!\n";
    return;
	}

	T *newarray = new T[goodSize];
	goodSize=0;
    minn = maxx = array[0];
    for(size_t i=0;i<size;i++) {
        if(mask[i]) {
            newarray[goodSize++] = array[i];
            if(array[i]<minn) minn=array[i];
            if(array[i]>maxx) maxx=array[i];
        }
    }

	mean = findMean(newarray,goodSize);
	stddev = findStddev(newarray,goodSize);
	median = findMedian(newarray,goodSize,true);
	madfm = findMADFM(newarray,goodSize,median,true);

	delete [] newarray;

}
template void findAllStats<short>(short*, size_t, bool*, short&, short&, short&,short&, short&, short&);
template void findAllStats<int>(int*, size_t, bool*, int&, int&, int &, int &, int &, int &);
template void findAllStats<long>(long*, size_t, bool*, long&, long&, long&, long&, long&, long&);
template void findAllStats<float>(float*, size_t, bool*, float&, float&, float&, float&, float&, float&);
template void findAllStats<double>(double*, size_t, bool*, double&, double&, double&, double&, double&, double&);



template <class T> 
T findMean(T *array, size_t size) {

  /// Find the mean of an array of numbers. Type independent.
  ///
  /// \param array 	The array of numbers.
  /// \param size 	The length of the array.
  /// \return 		The	mean value of the array, returned as a float.
  
	T mean = 0; 
	for(size_t i=0;i<size;i++) mean += array[i];
	if(size>0)
		mean /= T(size);
	return mean;
}
template short findMean<short>(short *, size_t);
template int findMean<int>(int *, size_t);
template long findMean<long>(long *, size_t);
template float findMean<float>(float *, size_t);
template double findMean<double>(double *, size_t);



template <class T>
T findMean(T *array, bool *mask, size_t size) {

  /// Find the mean of an array of numbers. Type independent.
  ///
  /// \param array	The array of numbers.
  /// \param mask 	An array of the same length that says whether to
  /// 		 		include each member of the array in the calculations. 
  ///		 		Only use values where mask=true. 
  /// \param size 	The length of the array.
  /// \return 		The mean value of the array, returned as a float.
  
	T mean = 0.;
	int ct=0;
	for(size_t i=0;i<size;i++){
		if(mask[i]){
			mean += array[i];
			ct++;
		}
	}
	if(ct>0) mean /= T(ct);
	return mean;
}
template short findMean<short>(short*, bool *, size_t);
template int findMean<int>(int *, bool *, size_t);
template long findMean<long>(long *, bool *, size_t);
template float findMean<float>(float *, bool *, size_t);
template double findMean<double>(double *, bool *, size_t);



template <class T> 
T findStddev(T *array, size_t size) {
 
  /// Find the rms or standard deviation of an array of
  /// numbers. Type independent. Calculated by iterating only once,
  /// using \sum x and \sum x^2 (no call to findMean).
  ///
  /// \param array 	The array of numbers.
  /// \param size 	The length of the array.
  /// \return 		The rms value of the array, returned as a float.

	T sumx=0,sumxx=0;
	T stddev=0;
	for(size_t i=0;i<size;i++){
		sumx += array[i];
		sumxx += (array[i]*array[i]);
	}
	if(size>0)
		stddev = sqrt(sumxx/T(size)-(sumx*sumx)/T(size*size));
	return stddev;
}
template short findStddev<short>(short*, size_t);
template int findStddev<int>(int *, size_t);
template long findStddev<long>(long *, size_t);
template float findStddev<float>(float *, size_t);
template double findStddev<double>(double *, size_t);



template <class T> 
T findStddev(T *array, bool *mask, size_t size) {
	
  /// Find the rms or standard deviation of an array of
  /// numbers. Type independent. Calculated by iterating only once,
  /// using \sum x and \sum x^2.
  ///
  /// \param array 	The array of numbers.
  /// \param mask 	An array of the same length that says whether to
  /// 		 		include each member of the array in the calculations. 
  ///		 		Only use values where mask=true.
  /// \param size 	The length of the array.
  /// \return 		The rms value of the array, returned as a float.

	T sumx=0,sumxx=0;
	T stddev=0;
	int ct=0;
	for(size_t i=0;i<size;i++){
		if(mask[i]){
		sumx += array[i];
		sumxx += (array[i]*array[i]);
		ct++;
		}
	}
	if(ct>0)
		stddev = sqrt(sumxx/T(ct)-(sumx*sumx)/T(ct*ct));
	return stddev;
}
template short findStddev<short>(short *, bool *, size_t);
template int findStddev<int>(int *, bool *, size_t);
template long findStddev<long>(long *, bool *, size_t);
template float findStddev<float>(float *, bool *, size_t);
template double findStddev<double>(double *, bool *, size_t);



template <class T>
T findMedian(T *array, size_t size, bool changeArray) {

  /// Find the median value of an array of numbers. Type independent.
  ///
  /// \param array 		 The array of numbers.
  /// \param size 		 The length of the array.
  /// \param changeArray Whether to use the provided array in calculations. 
  ///					 If true, the input array will be altered (ie. the order of 
  ///					 elements will be changed). Default value is "false".
  /// \return 			 The median value of the array, returned as the same type as 
  ///					 the array.
  
	T *newarray;
	if(changeArray) newarray = array;
	else{
		newarray = new T[size];
		for(size_t i=0;i<size;i++) newarray[i] = array[i];
	}
	T median;
	bool isEven = ((size%2)==0);
	std::nth_element(newarray,newarray+size/2,newarray+size);
	median = newarray[size/2];
	if(isEven){
		std::nth_element(newarray,newarray+size/2-1,newarray+size);
		median += newarray[size/2-1];
		median /= T(2);
	}
	if(!changeArray) delete [] newarray;
	return median;
}
template short findMedian<short>(short *, size_t, bool);
template int findMedian<int>(int *, size_t, bool);
template long findMedian<long>(long *, size_t, bool);
template float findMedian<float>(float *, size_t, bool);
template double findMedian<double>(double *, size_t, bool);


template <class T> 
T findMedian(std::vector<T> array, size_t size) {
	
	/// Find the median value of an array of numbers. Type independent.
  ///
  /// \param array 		 The array of numbers.
  /// \param size 		 The length of the array.
  /// \param changeArray Whether to use the provided array in calculations. 
  ///					 If true, the input array will be altered (ie. the order of 
  ///					 elements will be changed). Default value is "false".
  /// \return 			 The median value of the array, returned as the same type as 
  ///					 the array.
  
	T *newarray = new T[size];
	for(size_t i=0;i<size;i++) newarray[i] = array[i];
	
	T median;
	bool isEven = ((size%2)==0);
	std::nth_element(newarray,newarray+size/2,newarray+size);
	median = newarray[size/2];
	if(isEven){
		std::nth_element(newarray,newarray+size/2-1,newarray+size);
		median += newarray[size/2-1];
		median /= T(2);
	}
	delete [] newarray;
	return median;
}
template short findMedian<short>(std::vector<short>, size_t);
template int findMedian<int>(std::vector<int>, size_t);
template long findMedian<long>(std::vector<long>, size_t);
template float findMedian<float>(std::vector<float>, size_t);
template double findMedian<double>(std::vector<double>, size_t);


template <class T> 
T findMedian(T *array, bool *mask, size_t size) {
  
  /// Find the median value of an array of numbers. Type independent.
  ///
  /// \param array 	The array of numbers.
  /// \param size 	The length of the array.
  /// \param mask 	An array of the same length that says whether to
  /// 		 		include each member of the array in the calculations. 
  ///		 		Only use values where mask=true.
  /// \return 		The median value of the array, returned as the same type 
  ///				as the array.

	int goodSize=0,ct=0;
	for(size_t i=0;i<size;i++) if(mask[i]) goodSize++;
	T *newarray = new T[goodSize];
	for(size_t i=0;i<size;i++) {
		if(mask[i]) newarray[ct++] = array[i];
	}
	T median;
	bool isEven = ((goodSize%2)==0);
	std::nth_element(newarray,newarray+goodSize/2,newarray+goodSize);
	median = newarray[goodSize/2];
	if(isEven){
		std::nth_element(newarray,newarray+goodSize/2-1,newarray+goodSize);
		median += newarray[goodSize/2-1];
		median /= T(2);
	}
	delete [] newarray;
	return median;
}
template short findMedian<short>(short *, bool *, size_t);
template int findMedian<int>(int *, bool *, size_t);
template long findMedian<long>(long *, bool *, size_t);
template float findMedian<float>(float *, bool *, size_t);
template double findMedian<double>(double *, bool *, size_t);



template <class T> 
T findMADFM(T *array, size_t size, T median, bool changeArray) {

  /// Find the median absolute deviation from the median value of an
  /// array of numbers. Type independent. This version accepts a previously-
  /// calculated median value.
  /// 
  /// \param array 		 The array of numbers.
  /// \param size 		 The length of the array.
  /// \param median 	 The median of the array.
  /// \param changeArray Whether to use the provided array in calculations. 
  ///					 If true, the input array will be altered (ie. the order 
  ///					 of elements will be changed). Default value is "false".
  /// \return 			 The median absolute deviation from the median value of
  /// 		  			 the array, returned as the same type as the array.
  
	T *newarray;
	if(changeArray) newarray = array;
	else newarray = new T[size];

	T madfm;
	bool isEven = ((size%2)==0);
	for(size_t i=0;i<size;i++) newarray[i] = absval(array[i]-median);
	std::nth_element(newarray,newarray+size/2,newarray+size);
	madfm = newarray[size/2];
	if(isEven){
		std::nth_element(newarray,newarray+size/2-1,newarray+size);
		madfm += newarray[size/2-1];
		madfm /= T(2);
	}
	if(!changeArray) delete [] newarray;
	return madfm;
}
template short findMADFM<short>(short *, size_t, short, bool );
template int findMADFM<int>(int *, size_t, int, bool );
template long findMADFM<long>(long *, size_t, long, bool);
template float findMADFM<float>(float *, size_t, float, bool);
template double findMADFM<double>(double *, size_t, double, bool);



template <class T> 
T findMADFM(T *array, bool *mask, size_t size, T median) {

  /// Find the median absolute deviation from the median value of an
  /// array of numbers. Type independent. This version accepts a previously-
  /// calculated median value.
  /// 
  /// \param array 	The array of numbers.
  /// \param size 	The length of the array.
  /// \param median The median of the array.
  /// \param mask 	An array of the same length that says whether to
  /// 		 		include each member of the array in the calculations. 
  ///		 		Only use values where mask=true.
  /// \return 		The median absolute deviation from the median value of
  /// 		  		the array, returned as the same type as the array.
  
	int goodSize=0,ct=0;
	for(size_t i=0;i<size;i++) if(mask[i]) goodSize++;
	T *newarray = new T[goodSize];
	for(size_t i=0;i<size;i++){
		if(mask[i]) newarray[ct++] = absval(array[i]-median);
	}
	T madfm;
	bool isEven = ((goodSize%2)==0);
	std::nth_element(newarray,newarray+goodSize/2,newarray+goodSize);
	madfm = newarray[goodSize/2];
	if(isEven){
		std::nth_element(newarray,newarray+goodSize/2-1,newarray+goodSize);
		madfm += newarray[goodSize/2-1];
		madfm /= T(2);
	}
	delete [] newarray;
	return madfm;
}
template short findMADFM<short>(short*, bool *, size_t, short);
template int findMADFM<int>(int *, bool *, size_t, int);
template long findMADFM<long>(long *, bool *, size_t, long);
template float findMADFM<float>(float *, bool *, size_t, float);
template double findMADFM<double>(double *, bool *mask, size_t, double);
