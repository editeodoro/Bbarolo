//--------------------------------------------------------------------
// galfit_errors.cpp: Members functions of the Galfit class.
//--------------------------------------------------------------------

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
 along with BBarolo; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

 Correspondence concerning BBarolo may be directed to:
    Internet email: enrico.diteodoro@unibo.it
-----------------------------------------------------------------------*/


#include <iostream>
//#include <random>
#include "galfit.hh"
#include "galmod.hh"
#include "../Arrays/cube.hh"
#include "utils.hh"
#include "lsqfit.hh"
#include "progressbar.hh"

#define VROT  0
#define VDISP 1
#define DENS  2
#define Z0    3
#define INC   4
#define PA    5
#define XPOS  6
#define YPOS  7
#define VSYS  8
#define MAXPAR 9

using namespace std;
namespace Model {

template <class T>
void Galfit<T>::getErrors (Rings<T> *dr, T **err, int ir, T minimum) {
	
	for (int x=2; x--;) for (int y=nfree; y--;) err[x][y]=0.;
	
	ProgressBar bar(" Estimating errors... ", true);
    bar.setShowbar(in->pars().getShowbar());
	
	for (int ii=0; ii<2; ii++) {
		dr->vrot[ii]=outr->vrot[ir];
		dr->vdisp[ii]=outr->vdisp[ir];
		dr->dens[ii]=outr->dens[ir];
		dr->z0[ii]=outr->z0[ir];
		dr->inc[ii]=outr->inc[ir];
		dr->phi[ii]=outr->phi[ir];
		dr->xpos[ii]=outr->xpos[ir];
		dr->ypos[ii]=outr->ypos[ir];
		dr->vsys[ii]=outr->vsys[ir];
	}
			
	int free_var[nfree];
	float maxval[nfree],minval[nfree],midval[nfree];
	for (int nm=0, k=0; nm<MAXPAR; nm++) {
		if (mpar[nm]) {
			free_var[k]=nm;
			maxval[k]= -1.E04;
			minval[k]= 1.E04;
			if(nm==VROT) 		midval[k]=outr->vrot[ir];
			else if(nm==VDISP) 	midval[k]=outr->vdisp[ir];
			else if(nm==DENS) 	midval[k]=outr->dens[ir];
			else if(nm==Z0) 	midval[k]=outr->z0[ir];
			else if(nm==INC) 	midval[k]=outr->inc[ir];
			else if(nm==PA) 	midval[k]=outr->phi[ir];
			else if(nm==XPOS) 	midval[k]=outr->xpos[ir];
			else if(nm==YPOS) 	midval[k]=outr->ypos[ir];
			else if(nm==VSYS) 	midval[k]=outr->vsys[ir];
			k++;
		}
	}
	
	
	
	/*  METODO SERIO: TUTTI I PARAMETRI VARIANO CONTEMPORANEAMENTE  <<<---------------------------------
	int n_models=500*nfree;

	T *var_func = new T[n_models];
	T **var_val = allocate_2D<T>(n_models,nfree);	
			
	default_random_engine generator;	
	for (int nm=n_models;nm--;) {
		for (int nf=nfree; nf--;) {
			float minv, maxv;
			if (free_var[nf]==VROT) {
				//T maxvrot_r = fabs(AlltoVel<T>(in->getZphys(in->DimZ()-1),in->Head().Cunit(2), in->Head().Freq0())-dr->vsys[0]);
				//T maxvrot_l = fabs(AlltoVel<T>(in->getZphys(0),in->Head().Cunit(2),in->Head().Freq0())-dr->vsys[0]);
				maxv   = 200;//std::max(maxvrot_r*sin(dr->inc[0]*M_PI/180.),maxvrot_l*sin(dr->inc[0]*M_PI/180.));
				minv   = 0;//outr->vrot[ir]/2;
			}
			else {minv=mins[free_var[nf]];maxv=maxs[free_var[nf]];}
					

			//var_val[nm][nf] = unifrand(maxs[free_var[nf]], mins[free_var[nf]]);
			float stddev = min(maxv-midval[nf],midval[nf]-minv)/(1.386*2*sqrt(2));
			//normal_distribution<double> distribution(midval[nf],stddev);
			uniform_real_distribution<double> distribution(minv,maxv);
			var_val[nm][nf] = distribution(generator);
			//cout << midval[nf] << "  " << stddev <<  " " << var_val[nm][nf] << endl;

			switch(free_var[nf]) {
				case VROT:
					dr->vrot[0]=var_val[nm][nf]; 
					dr->vrot[1]=var_val[nm][nf];
					var_val[nm][nf]-=outr->vrot[ir];
					break;
				case VDISP:
					dr->vdisp[0]=var_val[nm][nf];
					dr->vdisp[1]=var_val[nm][nf];
					var_val[nm][nf]-=outr->vdisp[ir];
					break;
				case DENS:
					dr->dens[0]=var_val[nm][nf];
					dr->dens[1]=var_val[nm][nf];
					var_val[nm][nf]-=outr->dens[ir];
					break;
				case Z0:
					dr->z0[0]=var_val[nm][nf];
					dr->z0[1]=var_val[nm][nf];
					var_val[nm][nf]-=outr->z0[ir];
					break;
				case INC:
					dr->inc[0]=var_val[nm][nf];
					dr->inc[1]=var_val[nm][nf];
					var_val[nm][nf]-=outr->inc[ir];
					break;
				case PA:
					dr->phi[0]=var_val[nm][nf];
					dr->phi[1]=var_val[nm][nf];
					var_val[nm][nf]-=outr->phi[ir];
					break;
				case XPOS:
					dr->xpos[0]=var_val[nm][nf];
					dr->xpos[1]=var_val[nm][nf];
					var_val[nm][nf]-=outr->xpos[ir];
					break;
				case YPOS:
					dr->ypos[0]=var_val[nm][nf];
					dr->ypos[1]=var_val[nm][nf];
					var_val[nm][nf]-=outr->ypos[ir];
					break;
				case VSYS:
					dr->vsys[0]=var_val[nm][nf];
					dr->vsys[1]=var_val[nm][nf];
					var_val[nm][nf]-=outr->vsys[ir];
					break;
			}
					
			if (var_val[nm][nf]>maxval[nf]) maxval[nf] = var_val[nm][nf];
			if (var_val[nm][nf]<minval[nf]) minval[nf] = var_val[nm][nf];
										
		}
								
		var_func[nm] = model(dr);
		
		if (var_func[nm]<minimum) cout << dr->vrot[0] << "  " << dr->vdisp[0] << "  " << dr->inc[0] 
									   << "  " << dr->phi[0] << "   " << (var_func[nm]-minimum)/minimum << endl;
	}
			
			
			
	int n_bin = 40;
	for (int nf=nfree; nf--;) {
		std::string errname=in->pars().getOutfolder()+"errors"+to_string(ir)+"_"+to_string(free_var[nf])+".txt";
		std::ofstream errout(errname.c_str());
		float xrange = maxval[nf]-minval[nf];
		float bin_range =  xrange/n_bin;
		float xcounter = minval[nf];
		
		while (xcounter<=maxval[nf]) {
			float bin_value=10000;
			int nmod_inbin=0;
			for (int nm=n_models;nm--;) {
				if (var_val[nm][nf]<xcounter+bin_range &&
					var_val[nm][nf]>=xcounter) {
					float bin_val_temp = 100*(var_func[nm]-minimum)/minimum;
					if (bin_val_temp<bin_value) bin_value = bin_val_temp;
					//bin_value+=bin_val_temp;
					nmod_inbin++;
					}
						
				}
				//bin_value/=nmod_inbin;
				if (nmod_inbin==0) bin_value=0;
				errout << xcounter+bin_range/2. <<  "          " << bin_value << "  " << nmod_inbin << endl;
				xcounter+=bin_range;
			}
		errout.close();
	}
			
			
	deallocate_2D<T>(var_val,n_models);
	delete [] var_func;
	*/
			
			
	///*/// METODO ALL'ACQUA DI ROSE. UN PARAMETRO PER VOLTA
	
	
	int n_models=100;
	int n_bin = 40;
	if (verb) bar.init(n_models*nfree);
	uint cc=1;
	
	//default_random_engine generator;	
	for (int nf=nfree; nf--;) {
		for (int ii=0; ii<2; ii++) {
			dr->vrot[ii]=outr->vrot[ir];
			dr->vdisp[ii]=outr->vdisp[ir];
			dr->dens[ii]=outr->dens[ir];
			dr->z0[ii]=outr->z0[ir];
			dr->inc[ii]=outr->inc[ir];
			dr->phi[ii]=outr->phi[ir];
			dr->xpos[ii]=outr->xpos[ir];
			dr->ypos[ii]=outr->ypos[ir];
			dr->vsys[ii]=outr->vsys[ir];
		}

		T *var_func = new T[n_models];
		T *var_val  = new T[n_models];
		
		float minv, maxv;
		if (free_var[nf]==VROT) {
			//T maxvrot_r = fabs(AlltoVel<T>(in->getZphys(in->DimZ()-1),in->Head().Cunit(2), in->Head().Freq0())-dr->vsys[0]);
			//T maxvrot_l = fabs(AlltoVel<T>(in->getZphys(0),in->Head().Cunit(2),in->Head().Freq0())-dr->vsys[0]);
			//maxv   = 2*std::max(maxvrot_r*sin(dr->inc[0]*M_PI/180.),maxvrot_l*sin(dr->inc[0]*M_PI/180.));
			//minv   = outr->vrot[ir]/2;
			maxv = outr->vrot[ir]+outr->vrot[ir]/3.;
			minv = outr->vrot[ir]-outr->vrot[ir]/3.;
			if (minv<0) minv =0;
		}
		else if (free_var[nf]==VDISP) {
			maxv = outr->vdisp[ir]+15;
			minv = outr->vdisp[ir]-15;
			if (minv<0) minv =0;
		}		
		else {minv=mins[free_var[nf]];maxv=maxs[free_var[nf]];}
		 

		for (int nm=n_models;nm--;) {
			if (verb) bar.update(cc++);
			var_val[nm] = unifrand(maxv, minv);
			float stddev = 0.25509*min(maxv-midval[nf],midval[nf]-minv);    /// Read: 0.25509=1/(1.386*2*sqrt(2));
			//normal_distribution<double> distribution(midval[nf],stddev);
			//uniform_real_distribution<double> distribution(minv,maxv);
			//var_val[nm] = distribution(generator);
			//cout << midval[nf] << "  " << stddev <<  " " << var_val[nm] << endl;
					
			switch(free_var[nf]) {
				case VROT:
					dr->vrot[0]=var_val[nm]; 
					dr->vrot[1]=var_val[nm];
					var_val[nm]-=outr->vrot[ir];
					break;
				case VDISP:
					dr->vdisp[0]=var_val[nm];
					dr->vdisp[1]=var_val[nm];
					var_val[nm]-=outr->vdisp[ir];
					break;
				case DENS:
					dr->dens[0]=var_val[nm];
					dr->dens[1]=var_val[nm];
					var_val[nm]-=outr->dens[ir];
					break;
				case Z0:
					dr->z0[0]=var_val[nm];
					dr->z0[1]=var_val[nm];
					var_val[nm]-=outr->z0[ir];
					break;
				case INC:
					dr->inc[0]=var_val[nm];
					dr->inc[1]=var_val[nm];
					var_val[nm]-=outr->inc[ir];
					break;
				case PA:
					dr->phi[0]=var_val[nm];
					dr->phi[1]=var_val[nm];
					var_val[nm]-=outr->phi[ir];
					break;
				case XPOS:
					dr->xpos[0]=var_val[nm];
					dr->xpos[1]=var_val[nm];
					var_val[nm]-=outr->xpos[ir];
					break;
				case YPOS:
					dr->ypos[0]=var_val[nm];
					dr->ypos[1]=var_val[nm];
					var_val[nm]-=outr->ypos[ir];
					break;
				case VSYS:
					dr->vsys[0]=var_val[nm];
					dr->vsys[1]=var_val[nm];
					var_val[nm]-=outr->vsys[ir];
					break;
			}
					
			if (var_val[nm]>maxval[nf]) maxval[nf] = var_val[nm];
			if (var_val[nm]<minval[nf]) minval[nf] = var_val[nm];
					
			var_func[nm] = model(dr);
				
			//if (var_func[nm]<minimum) cout << dr->vrot[0] << "  " << dr->vdisp[0] << "  " << dr->inc[0] 
				//					   << "  " << dr->phi[0] << "   " << (var_func[nm]-minimum)/minimum << endl;
				
		}
		
		vector<T> xx_bin, yy_bin;
		float bin_range = (maxval[nf]-minval[nf])/n_bin;	
		//string errname=in->pars().getOutfolder()+"errors"+to_string(ir)+"_"+to_string(free_var[nf])+".txt";
		//ofstream errout(errname.c_str());
		for (int ib=0; ib<n_bin; ib++) {
			float bin_begin = minval[nf]+ib*bin_range;
			float bin_end   = bin_begin+bin_range;
			int nmod_inbin=0;
			xx_bin.push_back(bin_begin+bin_range/2.);
			yy_bin.push_back(10000);
			for (int nm=0; nm<n_models; nm++) {
				if (var_val[nm]>=bin_begin && var_val[nm]<bin_end) {
					T bin_val_temp = 100*(var_func[nm]-minimum)/minimum;
					if (bin_val_temp<yy_bin[ib]) yy_bin[ib] = bin_val_temp;
					//yy_bin[ib]+=bin_val_temp;
					nmod_inbin++;
				}
			}
			//yy_bin[ib]/=nmod_inbin;
			if (nmod_inbin==0) yy_bin[ib]=0;
			//errout << xx_bin[ib] <<  "    " << yy_bin[ib] << "    " << nmod_inbin << endl;
		}
		//errout.close();
		
		/// Fit with a polynomial with given order. Linear & costant term are set to 0 to impose the
		/// minimum in the best-fit model. If you want to change, you must change also after if (nrt<0) continue;
		/// where the roots for a cubic equation are calculated. 
		int order=3;
		T ww[n_bin], coeff[order+1], coefferr[order+1];
		bool mp[order+1];
		
		int n_datapoints=0, y_interv=15;
		while (n_datapoints<2*order) {
			for (int ib=n_bin; ib--;) {
				ww[ib] = (yy_bin[ib]>y_interv || yy_bin[ib]<=0 ) ? 0 : 1;
				n_datapoints+=ww[ib];
			}
			y_interv+=5;
		}
		
		//for (int ib=xx_bin.size(); ib--;) if (ww[ib]==1) errout << xx_bin[ib] <<  "    " << yy_bin[ib] << "    "  << endl;
		
		
		for (int i=order+1; i--;) {
			coeff[i]=coefferr[i]=0;
			mp[i] = true;
		}
		mp[0]=mp[1]=false;
		
		Lsqfit<T> lsq(&xx_bin[0],1,&yy_bin[0],ww,xx_bin.size(),coeff,coefferr,mp,order+1,&polyn,&polynd);
		int nrt = lsq.fit();

		//cout << setprecision(20) <<coeff[3] << "  " << coeff[2] << endl;
		
		
		xx_bin.clear();
		yy_bin.clear();	
		delete [] var_val;
		delete [] var_func;
		
		if (nrt<0) continue;
		
		///<<< We calculate errors within 5% of best model
		double threshold = 5.;   
		double a = coeff[2]/coeff[3];
		double b = coeff[1]/coeff[3];
		double c = -threshold/coeff[3];
		double Q = (a*a-3*b)/9.;
		double R = (2*a*a*a-9*a*b+27*c)/54.;
		
		double low_err=0, upp_err=0;
		if (R*R<Q*Q*Q) {
			// Three real roots
			double TT = acos(R/sqrt(Q*Q*Q));
			double x1 = -2*sqrt(Q)*cos(TT/3.)-a/3.;
			double x2 = -2*sqrt(Q)*cos((TT+2*M_PI)/3.)-a/3.;
			double x3 = -2*sqrt(Q)*cos((TT-2*M_PI)/3.)-a/3.;
			bool x1inb = x1>=minval[nf] && x1<=maxval[nf];
			bool x2inb = x2>=minval[nf] && x2<=maxval[nf];
			bool x3inb = x3>=minval[nf] && x3<=maxval[nf];
			bool allinb = x1inb && x2inb && x3inb;
			bool twoinb = !allinb && ((x1inb&&x2inb)||(x1inb&&x3inb)||(x2inb&&x3inb));
			bool oneinb = !allinb && !twoinb && (x1inb || x2inb || x3inb);
						
			if (allinb) {
				vector<double> rs(3);
				rs[0]=x1; rs[1]=x2; rs[2]=x3;
				sort(rs.begin(),rs.end());
				bool allPos = rs[0]>0 && rs[1]>0 && rs[2]>0; 
				bool twoPos = !allPos && ((rs[0]>0&&rs[1]>0)||(rs[0]>0&&rs[2]>0)||(rs[2]>0&&rs[1]>0)); 
				bool onePos = !allPos && !twoPos && (rs[0]>0 || rs[1]>0 || rs[2]>0);
				if (allPos) {low_err=-rs[0]; upp_err=rs[0];}
				else if (twoPos) {low_err=rs[0]; upp_err=rs[1];}
				else if (onePos) {low_err=rs[1]; upp_err=rs[2];}
				else {low_err=rs[2]; upp_err=-rs[2];}	
			}
			else if (twoinb) {
				double xx1=0, xx2=0;
				if (x1inb&&x2inb) {xx1=x1; xx2=x2;}
				else if (x1inb&&x3inb) {xx1=x1; xx2=x3;}
				else if (x2inb&&x3inb) {xx1=x2; xx2=x3;}
				if (xx1>xx2) swap(xx1,xx2);
				bool twoPos = xx1>0 && xx2>0;
				bool onePos = !twoPos && (xx1>0 || xx2>0);
				if (twoPos) {low_err=-xx1; upp_err=xx1;}
				else if (onePos) {low_err=xx1; upp_err=xx2;}
				else if (onePos) {low_err=xx2; upp_err=-xx2;}
			}
			else if (oneinb) {
				double xx1=0;
				if (x1inb) xx1=x1;
				else if (x2inb) xx1=x2;
				else if (x3inb) xx1=x3;
				if (xx1>0) {low_err=-xx1; upp_err=xx1;}
				else {low_err=xx1; upp_err=-xx1;}
			}
			else continue;
			
		}
		else {
			double A  = -cbrt(R+sqrt(R*R-Q*Q*Q));
			double B  = A!=0 ? Q/A : 0;			
			double x1 = (A+B)-a/3.;
			bool x1inb = x1>=minval[nf] && x1<=maxval[nf];
			if (!x1inb) continue;
			if (x1>0) {low_err=-x1; upp_err=x1;}
			else {low_err=x1; upp_err=-x1;}
		}
		
		err[0][nf] = low_err;
		err[1][nf] = upp_err;
	}
	//*/	
	if (verb) bar.fillSpace("Done.\n");

}
template void Galfit<float>::getErrors (Rings<float>*,float**,int,float);
template void Galfit<double>::getErrors (Rings<double>*,double**,int,double);

}

#undef VROT
#undef VDISP
#undef DENS
#undef Z0
#undef INC
#undef PA
#undef XPOS
#undef YPOS
#undef VSYS
#undef MAXPAR
