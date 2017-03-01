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

#ifndef SPACEPAR_HH_
#define SPACEPAR_HH_

#include <iostream> 
#include <cfloat> 
#include <Arrays/cube.hh>
#include <Utilities/galfit.hh>
#include <Utilities/utils.hh>
#include <Utilities/gnuplot.hh>
#include <Utilities/progressbar.hh>

using namespace Model;
template <class T>
class Spacepar : public Galfit<T>
{	
public:
	Spacepar(Cube<T> *c);
	~Spacepar(){};
	
	void calculate();

private:	

	std::string p1;
	std::string p2;
	T	minp1;
	T	minp2;
	T	maxp1;
	T	maxp2;
	T	delp1;
	T	delp2;
	
	void printmin(int n);
};


template <class T>
Spacepar<T>::Spacepar(Cube<T> *c) : Galfit<T>::Galfit(c) {
	
    p1 = Galfit<T>::in->pars().getP1();
    p2 = Galfit<T>::in->pars().getP2();
    p1 = makeupper(p1);
    p2 = makeupper(p2);

    if (p1!="VROT" && p1!="VDISP" && p1!="Z0"   && p1!="INC" &&
        p1!="PA"   && p1!="XPOS"  && p1!="YPOS" && p1!="VSYS") {
		std::cout << "SPACEPAR error: Unknown parameter 1 "
				  << p1 << std::endl;
		std::terminate();
	}
    if (p2!="VROT" && p2!="VDISP" && p2!="Z0"   && p2!="INC" &&
        p2!="PA"   && p2!="XPOS"  && p2!="YPOS" && p2!="VSYS") {
        std::cout << "SPACEPAR error: Unknown parameter 2 "
                  << p2 << std::endl;
        std::terminate();
    }
	if (p1==p2) {
		std::cout << "SPACEPAR error: p1=p2!!" << std::endl;
		std::terminate();
	}
	
	for (int i=0; i<6; i++) Galfit<T>::mpar[i]=false;

    minp1 = lround(Galfit<T>::in->pars().getP1p(0));
    minp2 = lround(Galfit<T>::in->pars().getP2p(0));
    maxp1 = lround(Galfit<T>::in->pars().getP1p(1));
    maxp2 = lround(Galfit<T>::in->pars().getP2p(1));
    delp1 = lround(Galfit<T>::in->pars().getP1p(2));
    delp2 = lround(Galfit<T>::in->pars().getP2p(2));
	
	Galfit<T>::maxs[0] = 10000;
	Galfit<T>::mins[0] = 0;
	Galfit<T>::maxs[1] = 10000;
	Galfit<T>::mins[1] = 0;
	Galfit<T>::maxs[3] = 10000;
	Galfit<T>::mins[3] = 0;
	Galfit<T>::maxs[4] = 90;
	Galfit<T>::mins[4] = 0;
	Galfit<T>::maxs[5] = 360;
	Galfit<T>::mins[5] = 0;
	Galfit<T>::maxs[6] = Galfit<T>::in->DimX()-1;
	Galfit<T>::mins[6] = 0;
	Galfit<T>::maxs[7] = Galfit<T>::in->DimY()-1;
	Galfit<T>::mins[7] = 0;

    std::string cunit2 = Galfit<T>::in->Head().Cunit(2);
    Galfit<T>::maxs[8] = AlltoVel(Galfit<T>::in->getZphys(Galfit<T>::in->DimZ()-1), Galfit<T>::in->Head());
    Galfit<T>::mins[8] = AlltoVel(Galfit<T>::in->getZphys(0), Galfit<T>::in->Head());
	
	if (Galfit<T>::in->pars().getSM()) {
		if (!Galfit<T>::setCfield()) {
			std::cout << "GALFIT warning: can not set an appropriate convolution "
					  << "field. Turning off the convolution step.\n";
			Galfit<T>::in->pars().setSM(false);
		}
	}
}


template <class T>
void Spacepar<T>::calculate() {
	
	bool verb = Galfit<T>::in->pars().isVerbose();
	
	if (verb) {	
		Galfit<T>::in->pars().setVerbosity(false);
		cout << showpoint << fixed << setprecision(2) << endl ;
		cout << setfill('=') << setw(40) << right << " SPACEBAR " << setw(34) << " ";
		cout << setfill(' ') << endl << endl;
		cout << " Varying "<<p1<<" in ["<<minp1<<","<<maxp1<<"] with del = "<<delp1
			 << " and "<<p2<<" in ["<<minp2<<","<<maxp2<<"] with del = "<<delp2;
		cout << endl << endl;
	}
	
	std::string filename = "./output/ring.dat";
	std::ofstream file;
	
    int start_rad = Galfit<T>::in->pars().getStartRad()<Galfit<T>::inr->nr ? Galfit<T>::in->pars().getStartRad() : 0;
    for (int ir=start_rad; ir<Galfit<T>::inr->nr; ir++) {
		
		file.open(filename.c_str());
		if (verb) cout << "\n\nRING NUMBER "<< ir+1 << ": \n";	

		int totmod = (maxp1-minp1)/(delp1)*(maxp2-minp2)/(delp2);
		totmod += (maxp1-minp1)/delp1+(maxp2-minp2)/delp2;
		
		int p1min=0, p2min=0;
        T funmin=FLT_MAX;
			
		int iminp1 = lround(minp1);
		int iminp2 = lround(minp2);
		int imaxp1 = lround(maxp1);
		int imaxp2 = lround(maxp2);
		int idelp1 = lround(delp1);
		int idelp2 = lround(delp2);
		
		int count=0;
        ProgressBar bar(" Calculating models...",true);
        bar.setShowbar(Galfit<T>::in->pars().getShowbar());
        if(verb) bar.init(totmod);

		for (int ip1=iminp1; ip1<=imaxp1; ip1+=idelp1) {
            for (int ip2=iminp2; ip2<=imaxp2; ip2+=idelp2) {
                if (verb) bar.update(count+1);

				Rings<T> *dring = new Rings<T>;
				Rings<T> *inri = Galfit<T>::inr;
				
				dring->nr = 2;
                float width1=0, width2=0;
                if (ir==0) width1 = width2 = (inri->radii[1]-inri->radii[0])/2.;
                else if (ir==inri->nr-1) width1 = width2 = (inri->radii[ir]-inri->radii[ir-1])/2.;
                else {
                    width1 = (inri->radii[ir]-inri->radii[ir-1])/2.;
                    width2 = (inri->radii[ir+1]-inri->radii[ir])/2.;
                }

                dring->radii.push_back(max(double(inri->radii[ir]-width1),0.));
                dring->radii.push_back(max(double(inri->radii[ir]+width2),0.));

                for (int i=0; i<dring->nr; i++) {
					if (p1=="VROT") dring->vrot.push_back(ip1);
					else if (p2=="VROT") dring->vrot.push_back(ip2);
					else dring->vrot.push_back(inri->vrot[ir]);
					
                    if (p1=="VDISP") dring->vdisp.push_back(ip1);
                    else if (p2=="VDISP") dring->vdisp.push_back(ip2);
                    else dring->vdisp.push_back(inri->vdisp[ir]);
					
					if (p1=="Z0") dring->z0.push_back(ip1);
					else if (p2=="Z0") dring->z0.push_back(ip2);
					else dring->z0.push_back(inri->z0[ir]);
					
					if (p1=="INC") dring->inc.push_back(ip1);
					else if (p2=="INC") dring->inc.push_back(ip2);
					else dring->inc.push_back(inri->inc[ir]);
					
					if (p1=="PA") dring->phi.push_back(ip1);
					else if (p2=="PA") dring->phi.push_back(ip2);
					else dring->phi.push_back(inri->phi[ir]);
					
					if (p1=="XPOS") dring->xpos.push_back(ip1);
					else if (p2=="XPOS") dring->xpos.push_back(ip2);
					else dring->xpos.push_back(inri->xpos[ir]);
					
					if (p1=="YPOS") dring->ypos.push_back(ip1);
					else if (p2=="YPOS") dring->ypos.push_back(ip2);
					else dring->ypos.push_back(inri->ypos[ir]);
					
					if (p1=="VSYS") dring->vsys.push_back(ip1);
					else if (p2=="VSYS") dring->vsys.push_back(ip2);
					else dring->vsys.push_back(inri->vsys[ir]);
					
					dring->dens.push_back(inri->dens[ir]);
				}

				T par[9];
				par[0] = dring->vrot.back();
				par[1] = dring->vdisp.back();
				par[2] = dring->dens.back();
				par[3] = dring->z0.back();
				par[4] = dring->inc.back();
				par[5] = dring->phi.back();
				par[6] = dring->xpos.back();		
				par[7] = dring->ypos.back();		
				par[8] = dring->vsys.back();		
				
                T minfunc = Galfit<T>::func3D(dring, par);
				
				if (minfunc<funmin) {
					funmin = minfunc;
					p1min = ip1;
					p2min = ip2;
				}
				
				file << ip1 << "  " << ip2 << "  " << minfunc << endl;
				
				delete dring;
				count++;
			}
			file <<	endl;
		}
		file.close();
		printmin(ir);

        if (verb) {
            bar.fillSpace("OK.\n");
			cout <<" Minimum ("<<funmin<<") at: "<<p1<<"="<<p1min<<" ,"<<p2<<"="<<p2min<<endl;
        }
	
	}
	
    remove(filename.c_str());
}


template <class T>
void Spacepar<T>::printmin (int nr) {

#ifdef HAVE_GNUPLOT
	Gnuplot gp;
	gp.begin();	
    gp.commandln("set terminal png enhanced");
    //gp.commandln("set terminal postscript eps color enhanced");
    gp.commandln("set pm3d map");
	gp.commandln("set size square");
	gp.commandln("unset key");
	std::string cmd = "set xrange ["+to_string(minp1)+":"+to_string(maxp1)+"]";
	gp.commandln(cmd.c_str());
	cmd = "set yrange ["+to_string(minp2)+":"+to_string(maxp2)+"]";	
	gp.commandln(cmd.c_str());
    cmd = Galfit<T>::in->pars().getOutfolder()+"r"+to_string<int>(nr+1)+".png";
	gp.commandln(("set output '"+cmd+"'").c_str());
    if 		 (p1=="VROT") 	cmd="set xlabel 'V_{rot}  [Km/s]'";
	else if (p1=="VDISP") cmd="set xlabel 'Dispersion [km/s]'";
	else if (p1=="Z0") 	cmd="set xlabel 'Scale height [arcs]'";
	else if (p1=="INC") 	cmd="set xlabel 'Inclination [degree]'";
	else if (p1=="PA") 	cmd="set xlabel 'Position angle [degree]'";
	else if (p1=="XPOS") 	cmd="set xlabel 'Xcenter [pixel]'";
	else if (p1=="YPOS") 	cmd="set xlabel 'Ycenter [pixel]'";
	else if (p1=="VDISP") 	cmd="set xlabel 'Systemic velocity [km/s]'";
	gp.commandln(cmd.c_str());
	if 		 (p2=="VROT")  cmd="set ylabel 'V_rot  [Km/s]'";
	else if (p2=="VDISP") cmd="set ylabel 'Dispersion [km/s]'";
	else if (p2=="Z0") 	cmd="set ylabel 'Scale height [arcs]'";
	else if (p2=="INC") 	cmd="set ylabel 'Inclination [degree]'";
	else if (p2=="PA") 	cmd="set ylabel 'Position angle [degree]'";
	else if (p1=="XPOS") 	cmd="set xlabel 'Xcenter [pixel]'";
	else if (p1=="YPOS") 	cmd="set xlabel 'Ycenter [pixel]'";
	else if (p1=="VDISP") 	cmd="set xlabel 'Systemic velocity [km/s]'";
	gp.commandln(cmd.c_str());
	
	cmd = "set title 'Ring #"+to_string(nr+1)+" ("+
		  to_string(Galfit<T>::inr->radii[nr])+" arcs)'";
	gp.commandln(cmd.c_str());
	
	cmd = "splot './output/ring.dat'";
	gp.commandln(cmd.c_str());

	gp.end();
#endif	
}

#endif
