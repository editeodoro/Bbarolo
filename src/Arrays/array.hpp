//------------------------------------------------------------------
// array.hh
// Array class
//------------------------------------------------------------------  
#include <iostream>

#ifndef ARRAY_HH
#define ARRAY_HH

/////////////////////////////////////////////////////////////////////
// N Dimentional Array Class
/////////////////////////////////////////////////////////////////////
template<class Type, size_t N>
class Array {
public:

    // BASIC CONSTRUCTORS
    Array() {Ndim=N;}
    Array(size_t *Dim) {setArray(Dim);}
    Array(size_t *Dim, Type v) {setArray(Dim, v);}
    Array(size_t *Dim, Type *vp) {setArray(Dim,vp);}
    Array(const Array& a) {this->operator=(a);};
    ~Array() {if(p) delete[] p;}

    // SPECIAL CONSTRUCTORS FOR 4D,3D,2D and 1D ARRAYS
    Array(size_t xi,size_t yi, size_t zi, size_t qi) {if (checkNDim(4)) setArray(xi,yi,zi,qi);}
    Array(size_t xi,size_t yi, size_t zi) {if (checkNDim(3)) setArray(xi,yi,zi);}
    Array(size_t xi,size_t yi) {if (checkNDim(2)) setArray(xi,yi);}
    Array(size_t xi) {if (checkNDim(1)) setArray(xi);}

    void setArray(size_t *Dim) {
        Ndim=N; size=1;
        for (int i=N; i--;) {dim[i]=Dim[i]; size*=dim[i];}
        if (p) delete [] p;
        p=new Type[size];
    }

    void setArray(size_t *Dim, Type v) {setArray(Dim); for (int i=size; i--;) p[i]=v;}
    void setArray(size_t *Dim, Type *v) {setArray(Dim); for (int i=size; i--;) p[i]=v[i];}


    // OPERATOR OVERLOADING
    Array& operator=(const Array& a) {
        if (checkNDim(a.Ndim)) {
            size_t Dim[N];
            for (int i=N; i--;) Dim[i]=a.dim[i];
            setArray(Dim,a.p);
        }
        return *this;
    }

    inline Type& operator()(size_t *pix) {
#ifdef INDEXCHECK
        bool goodPix = true;
        for (int i=N; i--;) goodPix = goodPix && (pix[i]<dim[i]);
        if(!goodPix) {
            std::cerr << "\nBad index (";
            for (int i=N; i--;) std::cerr << pix[i] << ",";
            std::cerr << ") > (";
            for (int i=N; i--;) std::cerr << dim[i] << ",";
            std::cerr << ")\n";
            std::terminate();
        }
#endif

        size_t pixel = nPix(dim);
        return *(p+pixel);

    }

    inline Type& operator()(size_t i) {return *(p+i);}
    inline Type& operator[] (size_t s) {return this->operator()(s);}

    // FRIEND FUNCTIONS
    template <class T, size_t M, class K> friend Array<T,M>& operator+=(Array<T,M>& a,const K& b);
    template <class T, size_t M, class K> friend Array<T,M>& operator-=(Array<T,M>& a,const K& b);
    template <class T, size_t M, class K> friend Array<T,M>& operator*=(Array<T,M>& a,const K& b);
    template <class T, size_t M, class K> friend Array<T,M>& operator/=(Array<T,M>& a,const K& b);

    template <class T, size_t M> friend Array<T,M>& operator+=(Array<T,M>& a,const Array<T,M>& b);
    template <class T, size_t M> friend Array<T,M>& operator-=(Array<T,M>& a,const Array<T,M>& b);
    template <class T, size_t M> friend Array<T,M>& operator*=(Array<T,M>& a,const Array<T,M>& b);
    template <class T, size_t M> friend Array<T,M>& operator/=(Array<T,M>& a,const Array<T,M>& b);

    template <class T, size_t M> friend bool operator==(Array<T,M>& a, Array<T,M>& b);
    template <class T, size_t M> friend bool operator!=(Array<T,M>& a, Array<T,M>& b);


    // OTHER UTILITY FUNCTIONS
    inline size_t nPix (size_t *pix) {
        size_t pixel=0;
        for (int i=0; i<N; i++) {
            size_t dimf = 1;
            for (int n=0; n<i; n++) dimf*=dim[n];
            pixel += (pix[i]*dimf);
        }
        return pixel;
    }


    size_t Size () {return size;}
    size_t* Dim () {return dim;}
    size_t& Dim (size_t n) {return dim[n];}
    Type* P () {return p;}
    Type& P (size_t i) {return p[i];}


    // SPECIALIZED FUNCTION TO BE USED FOR 4D, 3D and 2D ARRAYS
    // WARNING: no check on pixels or dimensions ---> use consciously
    void setArray(size_t xi,size_t yi, size_t zi, size_t qi) {size_t Dim[4] = {xi,yi,zi,qi}; setArray(Dim);}
    void setArray(size_t xi,size_t yi, size_t zi) {size_t Dim[3] = {xi,yi,zi}; setArray(Dim);}
    void setArray(size_t xi,size_t yi) {size_t Dim[2] = {xi,yi}; setArray(Dim);}
    inline Type& operator()(size_t i,size_t j,size_t k, size_t l) {return *(p+i+j*dim[0]+k*dim[0]*dim[1]+l*dim[0]*dim[1]*dim[2]);}
    inline Type& operator()(size_t i,size_t j,size_t k) {return *(p+i+j*dim[0]+k*dim[0]*dim[1]);}
    inline Type& operator()(size_t i,size_t j) {return *(p+i+j*dim[0]);}
    inline size_t nPix (size_t i,size_t j,size_t k, size_t l) {return (i+j*dim[0]+k*dim[0]*dim[1]+l*dim[0]*dim[1]*dim[2]);}
    inline size_t nPix (size_t i,size_t j,size_t k) {return (i+j*dim[0]+k*dim[0]*dim[1]);}
    inline size_t nPix (size_t i,size_t j) {return (i+j*dim[0]);}


private:
    Type* 	p = NULL;
    size_t 	dim[N], Ndim, size;


    bool checkSize (const Array &a, bool terminate=true) const {
        if (size!=a.size) {
            if (terminate) {
                std::cerr << "\nError: arrays must have the same size: "
                          << size << " != " << a.size << std::endl;
                std::terminate();
            }
            else return false;
        }

        bool sameSize = true;
        for (int i=N; i--;) sameSize = sameSize && (dim[i]==a.dim[i]);
        if(!sameSize) {
            if (terminate) {
                std::cerr << "\nError: arrays must have the same shape: (";
                for (int i=N; i--;) std::cerr << dim[i] << ",";
                std::cerr << ") != (";
                for (int i=N; i--;) std::cerr << a.dim[i] << ",";
                std::cerr << ")\n";
                std::terminate();
            }
            else return false;
        }
        else return true;
    }

    bool checkNDim (size_t nd) {
        if (nd!=N) {
            std::cerr << "\nError: you must provide "<< N << " dimensions. \n";
            std::terminate();
        }
        else return true;
    }
};


// OPERATIONS WITH SCALARS
template <class T, size_t M, class K>
Array<T,M>& operator+=(Array<T,M>& a, const K& b) {for (int i=a.size; i--;) a.p[i] += T(b);}

template <class T, size_t M, class K>
Array<T,M>& operator-=(Array<T,M>& a, const K& b) {for (int i=a.size; i--;) a.p[i] -= T(b);}
template <class T, size_t M>
Array<T,M>& operator-=(Array<T,M>&& a, double& b);

template <class T, size_t M, class K>
Array<T,M>& operator*=(Array<T,M>& a, const K& b) {for (int i=a.size; i--;) a.p[i] *= T(b);}

template <class T, size_t M, class K>
Array<T,M>& operator/=(Array<T,M>& a, const K& b) {for (int i=a.size; i--;) a.p[i] /= T(b);}

template <class T, size_t M, class K>
Array<T,M> operator+(const Array<T,M>& a, const K& b) {Array<T,M> c=a; c+=b; return c;}

template <class T, size_t M, class K>
Array<T,M> operator-(const Array<T,M>& a, const K& b) {Array<T,M> c=a; c-=b; return c;}

template <class T, size_t M, class K>
Array<T,M> operator*(const Array<T,M>& a, const K& b) {Array<T,M> c=a; c*=b; return c;}

template <class T, size_t M, class K>
Array<T,M> operator/(const Array<T,M>& a, const K& b) {Array<T,M> c=a; c/=b; return c;}

template <class T, size_t M, class K>
Array<T,M> operator+(const K& b, const Array<T,M>& a) {return a+b;}

template <class T, size_t M, class K>
Array<T,M> operator-(const K& b, const Array<T,M>& a) {return a-b;}

template <class T, size_t M, class K>
Array<T,M> operator*(const K& b, const Array<T,M>& a) {return a*b;}

template <class T, size_t M, class K>
Array<T,M> operator/(const K& b, const Array<T,M>& a) {return a/b;}


// OPERATIONS WITH OTHER ARRAYS
template <class T, size_t M>
Array<T,M>& operator+=(Array<T,M>& a,const Array<T,M>& b) {
    if (a.checkSize(b)) for (int i=a.size; i--;) a.p[i] += b.p[i];
    return a;
}

template <class T, size_t M>
Array<T,M>& operator-=(Array<T,M>& a,const Array<T,M>& b) {
    if (a.checkSize(b)) for (int i=a.size; i--;) a.p[i] -= b.p[i];
    return a;
}

template <class T, size_t M>
Array<T,M>& operator*=(Array<T,M>& a,const Array<T,M>& b) {
    if (a.checkSize(b)) for (int i=a.size; i--;) a.p[i] *= b.p[i];
    return a;
}

template <class T, size_t M>
Array<T,M>& operator/=(Array<T,M>& a,const Array<T,M>& b) {
    if (a.checkSize(b)) for (int i=a.size; i--;) a.p[i] /= a.p[i];
    return a;
}	

template <class T, size_t M>
Array<T,M> operator+(const Array<T,M>& a, const Array<T,M>& b) {Array<T,M> c=a; return (c+=b);}

template <class T, size_t M>
Array<T,M> operator-(const Array<T,M>& a, const Array<T,M>& b) {Array<T,M> c=a; return (c+=b);}

template <class T, size_t M>
Array<T,M> operator*(const Array<T,M>& a, const Array<T,M>& b) {Array<T,M> c=a; return (c+=b);}

template <class T, size_t M>
Array<T,M> operator/(const Array<T,M>& a, const Array<T,M>& b) {Array<T,M> c=a; return (c+=b);}


// OTHER OPERATORS
template <class T, size_t M>
bool operator==(const Array<T,M>& a, const Array<T,M>& b) {
    if (!a.checkSize(b,false)) return false;
    for(int i=a.size; i--;) if(a.p[i]!=b.p[i]) return false;
    return true;
}

template <class T, size_t M>
bool operator!=(const Array<T,M>& a, const Array<T,M>& b) {return(!(a==b));}


typedef Array<double,4> double4D;
typedef Array<float,4>  float4D;
typedef Array<long,4>   long4D;
typedef Array<int,4>    int4D;
typedef Array<bool,4>   bool4D;
typedef Array<short,4>  short4D;
typedef Array<double,3> double3D;
typedef Array<float,3>  float3D;
typedef Array<long,3>   long3D;
typedef Array<int,3>    int3D;
typedef Array<bool,3>   bool3D;
typedef Array<short,3>  short3D;
typedef Array<double,2> double2D;
typedef Array<float,2>  float2D;
typedef Array<long,2>   long2D;
typedef Array<int,2>    int2D;
typedef Array<short,2>  short2D;
typedef Array<bool,2>   bool2D;


#endif
