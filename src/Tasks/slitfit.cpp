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

#include <Tasks/galfit.hh>
#include <Arrays/cube.hh>
#include <Arrays/image.hh>
#include <Utilities/utils.hh>
#include <Utilities/gnuplot.hh>



namespace Model {

template <class T>
void Galfit<T>::slit_init(Cube<T> *c) {

    Param &p = c->pars();
    par = p.getParGF();


    // Creating a dummy datacube needed in input for the fit
    int axisdim[3] = {c->DimY(),c->DimY(),c->DimX()};
    in = new Cube<T>(axisdim);
    in->saveParam(c->pars());

    Header &h = in->Head();
    h.setNumAx(3);
    h.setDimAx(0,axisdim[0]);
    h.setDimAx(1,axisdim[1]);
    h.setDimAx(2,axisdim[2]);

    h.setCrpix(0,in->DimX()/2.+1);
    h.setCrval(0,0.);
    h.setCdelt(0,-fabs(c->Head().Cdelt(1)));
    h.setCunit(0,c->Head().Cunit(1));
    h.setCtype(0,"RA---SIN");
    h.setCrpix(1,in->DimY()/2.+1);
    h.setCrval(1,0.);
    h.setCdelt(1,fabs(c->Head().Cdelt(1)));
    h.setCunit(1,c->Head().Cunit(1));
    h.setCtype(1,"DEC--SIN");
    h.setCrpix(2,c->Head().Crpix(0));
    h.setCrval(2,c->Head().Crval(0));
    h.setCdelt(2,c->Head().Cdelt(0));
    h.setCunit(2,c->Head().Cunit(0));
    h.setCtype(2,c->Head().Ctype(0));

    h.setName(c->Head().Name());

    h.setFreq0(0.1420405751786E10);
    h.setMinMax(0.,0.);
    in->saveHead(h);

    Rings<T> *inR = readRings<T>(c->pars().getParGF(),c->Head());

    // Setup all needed parameters
    bool isVerbose = p.isVerbose();
    in->pars().setVerbosity(false);
    setup(in,inR,&par);
    in->pars().setVerbosity(isVerbose);


    wpow = 0;
    par.NORMALCUBE = false;
    par.NV = 100;

    func_norm = &Model::Galfit<T>::slitfunc;

    showInitial(inr, std::cout);


    // Now defining the line images
    line_im = new Cube<T>(c->AxisDim());
    line_imDefined = true;
    line_im->saveHead(c->Head());
    line_im->saveParam(c->pars());
    for (int i=0; i<c->NumPix(); i++) line_im->Array(i) = c->Array(i);


/*
    // Continuum subctration.
    int w = 11;
    for (int y=0; y<line_im->DimY(); y++) {
        T val0 = findMean(&line_im->Array()[line_im->nPix(0,y,0)],2);
        T val1 = findMean(&line_im->Array()[line_im->nPix(line_im->DimX()-w-1,y,0)],w);
        if (val0==0 && val1==0) continue;
        int x0 = floor(w/2.);
        int x1 = line_im->DimX()-ceil(w/2.);
        T slope = (val1-val0)/(x1-x0);
        T intercept = -slope*x0+val0;
        for (int x=0;x<line_im->DimX();x++) (*line_im)(x,y,0)-= (slope*x+intercept);
    }
*/

    // Identification of the line
    line_im->pars().getParSE().minChannels = 1;
    line_im->pars().getParSE().minPix = 10;
    line_im->search();
    int numObj = line_im->getNumObj();
    if (numObj==0)  {
        std::cout << "SLITFIT error: No lines detected in the datacube. Cannot fit!!! \n";
        exit(EXIT_FAILURE);
    }
    else if (numObj>1) {
        uint n=0, size=0;
        for (int i=0; i<numObj; i++)
            if (line_im->pObject(i)->getSize()>size) {n=i;size=line_im->pObject(i)->getSize();}
        for (int i=0; i<numObj; i++)
            if (i!=n) line_im->getSources()->pObjectList()->erase(line_im->getSources()->pObjectList()->begin()+i);
    }
/*
    std::vector<bool> isObj(line_im->NumPix(),false);
    typename std::vector<Voxel<T> > voxelList = line_im->getSources()->LargestDetection()->getPixelSet(line_im->Array(), line_im->AxisDim());
    typename std::vector<Voxel<T> >::iterator vox;
    for(vox=voxelList.begin();vox<voxelList.end();vox++)
        isObj[line_im->nPix(vox->getX(),vox->getY(),vox->getZ())] = true;

    for (int i=0; i<line_im->NumPix(); i++) line_im->Array(i)*=isObj[i];
*/
    line_im->fitswrite_3d((p.getOutfolder()+in->Head().Name()+".fits").c_str());







      /*
    line_im->Head().setCrpix(0,h.Crpix(2));
    line_im->Head().setCrpix(1,h.Crpix(1));
    line_im->Head().setCrval(0,h.Crval(2));
    line_im->Head().setCrval(1,0.);
    line_im->Head().setCdelt(0,h.Cdelt(2));
    line_im->Head().setCdelt(1,h.Cdelt(1));
    line_im->Head().setCunit(0,h.Cunit(2));
    line_im->Head().setCunit(1,h.Cunit(1));
    line_im->Head().setCtype(0,h.Ctype(2));
    line_im->Head().setCtype(1,"Offset");
    line_im->Head().setMinMax(0.,0.);

    in = fcube;



    in->search();

    int numObj = in->getNumObj();
    if (numObj==0)  {
        std::cout << "SLITFIT error: No lines detected in the datacube. Cannot fit!!! \n";
        exit(EXIT_FAILURE);
    }
    else if (numObj>1) {
        uint n=0, size=0;
        for (int i=0; i<numObj; i++)
            if (in->pObject(i)->getSize()>size) {n=i;size=in->pObject(i)->getSize();}
        for (int i=0; i<numObj; i++)
            if (i!=n) in->getSources()->pObjectList()->erase(in->getSources()->pObjectList()->begin()+i);
    }
    
    
    // Setting wavelegths
    Image2D<float> *Wave = new Image2D<float>;
    Wave->Head().setWarning(false);
    Wave->readImage(in->pars().getWavefile());
    if (in->DimX()!=Wave->DimX()) {
        std::cout << "SLITFIT error: The data and the wave files have different dimensions. \n";
        exit(EXIT_FAILURE);
    }

    double *wave = new double[Wave->DimX()];
    for (int i=0; i<Wave->DimX(); i++) wave[i]=Wave->Array(i);
    delete Wave;

    double wave_s=0, wave_rest=0;
    if (p->getLine()=="Ha") wave_rest = Ha;
    else if (p->getLine()=="Hb") wave_s = Hb;
    else {
        std::cout << "SLITFIT: Sorry, I do not know the line " << p->getLine() << std::endl;
        exit(EXIT_FAILURE);
    }
    wave_s = wave_rest*(p->getRedshift()+1);

    if (wave_s<wave[0] || wave_s>wave[in->DimX()-1]) {
        std::cout << "SLITFIT error: the line is not inside the provided wavelength range." << std::endl;
        exit(EXIT_FAILURE);
    }

    int chan_s=0;
    double diff_min=1.E06;
    for (int i=0;i<in->DimX();i++) {
        double diff_w = fabs(wave[i]-wave_s);
        if (diff_w<diff_min) {diff_min=diff_w; chan_s=i;};
    }

    //std::cout << in->LargestDetection()->getYaverage()<<"  " << wave_s<< " "<< chan_s <<std::endl;


    // Setting spectral broadedning
    double checkIvar = strtod(in->pars().getIvarfile().c_str(),NULL);
    // TO DO STUFF WITH IVAR FILE


    // Set the image of the line and continuum subtraction.
    int offset = 50;
    int spat_s = lround(in->getSources()->LargestDetection()->getYaverage());
    int minx = (chan_s-offset)>=0 ? chan_s-offset : 0;
    int maxx = (chan_s+offset)<in->DimX() ? chan_s+offset : in->DimX()-1;
    int miny = (spat_s-offset)>=0 ? spat_s-offset : 0;
    int maxy = (spat_s+offset)<in->DimY() ? spat_s+offset : in->DimY()-1;

    int axisdim[3] = {maxx-minx+1,maxy-miny+1,1};
    line_im = new Cube<T>(axisdim);
    line_imDefined = true;
    line_im->saveHead(in->Head());
    line_im->saveParam(in->pars());
    for (int x=0; x<line_im->DimX();x++)
        for (int y=0; y<line_im->DimY();y++)
           (*line_im)(x,y,0) = (*in)(x+minx,y+miny,0);


    // Continuum subctration.
    int w = 11;
    for (int y=0; y<line_im->DimY(); y++) {
        T val0 = findMean(&line_im->Array()[line_im->nPix(0,y,0)],2);
        T val1 = findMean(&line_im->Array()[line_im->nPix(line_im->DimX()-w-1,y,0)],w);
        if (val0==0 && val1==0) continue;
        int x0 = floor(w/2.);
        int x1 = line_im->DimX()-ceil(w/2.);
        T slope = (val1-val0)/(x1-x0);
        T intercept = -slope*x0+val0;
        for (int x=0;x<line_im->DimX();x++) (*line_im)(x,y,0)-= (slope*x+intercept);
    }


    // Identification of the line
    line_im->search();
    numObj = line_im->getNumObj();
    if (numObj==0)  {
        std::cout << "SLITFIT error: No lines detected in the datacube. Cannot fit!!! \n";
        exit(EXIT_FAILURE);
    }
    else if (numObj>1) {
        uint n=0, size=0;
        for (int i=0; i<numObj; i++)
            if (line_im->pObject(i)->getSize()>size) {n=i;size=line_im->pObject(i)->getSize();}
        for (int i=0; i<numObj; i++)
            if (i!=n) line_im->getSources()->pObjectList()->erase(line_im->getSources()->pObjectList()->begin()+i);
    }

    std::vector<bool> isObj(line_im->NumPix(),false);
    typename std::vector<Voxel<T> > voxelList = line_im->getSources()->LargestDetection()->getPixelSet(line_im->Array(), line_im->AxisDim());
    typename std::vector<Voxel<T> >::iterator vox;
    for(vox=voxelList.begin();vox<voxelList.end();vox++)
        isObj[line_im->nPix(vox->getX(),vox->getY(),vox->getZ())] = true;

    for (int i=0; i<line_im->NumPix(); i++) (*line_im)(i)*=isObj[i];


    ///////////////////////////////////

    arcconv = arcsconv(in->Head().Cunit(0));
    distance = RedtoDist(p->getRedshift());

    if (in->Head().BeamArea()==0) {
        cout << "\n Beam information is not available in the header: assuming a "
             << in->pars().getBeamFWHM()*3600 << " arcsec beam. \n You can set the beam "
             << "with BeamFWHM parameter (in arcsec).\n\n";
        in->Head().setBmaj(in->pars().getBeamFWHM());
        in->Head().setBmin(in->pars().getBeamFWHM());
        in->Head().calcArea();
    }
    // Try to read ring information from an input file
    Rings<T> file_rings;
    bool radii_b,xpos_b,ypos_b,vsys_b,vrot_b,vdisp_b,z0_b,dens_b,inc_b,pa_b;
    radii_b = getDataColumn(file_rings.radii,par.RADII);
    xpos_b  = getDataColumn(file_rings.xpos,par.XPOS);
    ypos_b  = getDataColumn(file_rings.ypos,par.YPOS);
    vsys_b  = getDataColumn(file_rings.vsys,par.VSYS);
    vrot_b  = getDataColumn(file_rings.vrot,par.VROT);
    vdisp_b = getDataColumn(file_rings.vdisp,par.VDISP);
    z0_b    = getDataColumn(file_rings.z0,par.Z0);
    dens_b  = getDataColumn(file_rings.dens,par.DENS);
    inc_b   = getDataColumn(file_rings.inc,par.INC);
    pa_b    = getDataColumn(file_rings.phi,par.PHI);

    size_t size[10] = {file_rings.radii.size(),file_rings.xpos.size(),
                    file_rings.ypos.size(), file_rings.vsys.size(),
                    file_rings.vrot.size(),file_rings.vdisp.size(),
                    file_rings.z0.size(),file_rings.dens.size(),
                    file_rings.inc.size(),file_rings.phi.size()};

    int max_size=INT_MAX;
    for (int i=0; i<10; i++) if (size[i]!=0 && size[i]<max_size) max_size=size[i];

    int nr=0;
    T radsep, xpos, ypos, vsys, vrot, vdisp, z0, dens, inc, pa;
    nr    = par.NRADII;
    radsep= par.RADSEP;
    vrot  = atof(par.VROT.c_str());
    inc   = atof(par.INC.c_str());
    vdisp = par.VDISP!="-1" ? atof(par.VDISP.c_str()): 8.;                  // default is 8 km/s
    z0    = par.Z0!="-1" ? atof(par.Z0.c_str()) : 0.15/KpcPerArc(distance); // default is 150 parsec
    dens  = par.DENS!="-1" ? atof(par.DENS.c_str()) : 1.;
    xpos  = par.XPOS!="-1" ? atof(par.XPOS.c_str()) : in->Head().Crpix(0)-1;
    ypos  = par.YPOS!="-1" ? atof(par.YPOS.c_str()) : in->Head().Crpix(1)-1;
    vsys  = par.VSYS!="-1" ? atof(par.VSYS.c_str()) : 0.;
    pa    = par.PHI!="-1"  ? atof(par.PHI.c_str()) : 0;

    if (pa!=180 && pa!=0) {
        std::cout << "SLITFIT WARNING: PA must be 0 or 180. Setting to 0.\n";
        pa = 0;
    }

    nr = nr>0 && nr<max_size ? nr : max_size;
    if (radii_b) {
       radsep = 0;
       for (uint i=1; i<file_rings.radii.size()-1; i++)
            radsep += file_rings.radii[i+1]-file_rings.radii[i];
        radsep/=(file_rings.radii.size()-2);
     }

     inr = new Rings<T>;
     inDefined = true;
     inr->nr    = nr;
     inr->radsep = radsep;
     for (int i=0; i<inr->nr; i++) {
        if (radii_b) inr->radii.push_back(file_rings.radii[i]);
        else inr->radii.push_back(i*radsep+radsep/2.);
        if (vrot_b) inr->vrot.push_back(file_rings.vrot[i]);
        else inr->vrot.push_back(vrot);
        if (vdisp_b) inr->vdisp.push_back(file_rings.vdisp[i]);
        else inr->vdisp.push_back(vdisp);
        if (z0_b) inr->z0.push_back(file_rings.z0[i]);
        else inr->z0.push_back(z0);
        if (dens_b) inr->dens.push_back(file_rings.dens[i]*1.E20);
        else inr->dens.push_back(dens*1.E20);
        if (inc_b) inr->inc.push_back(file_rings.inc[i]);
        else inr->inc.push_back(inc);
        if (pa_b) inr->phi.push_back(file_rings.phi[i]);
        else inr->phi.push_back(pa);
        if (xpos_b) inr->xpos.push_back(file_rings.xpos[i]);
        else inr->xpos.push_back(xpos);
        if (ypos_b) inr->ypos.push_back(file_rings.ypos[i]);
        else inr->ypos.push_back(ypos);
        if (vsys_b) inr->vsys.push_back(file_rings.vsys[i]);
        else inr->vsys.push_back(vsys);
        inr->vrad.push_back(0);
        
        // In the current version, vertical motions, and gradients are not fitted
        inr->vvert.push_back(0);
        inr->dvdz.push_back(0);
        inr->zcyl.push_back(0);
        
    }

    if (inr->radii[0]!=0) inr->radii[0]=0;

    setFree();


    setup(in, inr, &par);

    outr = new Rings<T>;
    *outr = *inr;
    outDefined = true;

    func_norm = &Model::Galfit<T>::slitfunc;

    delete [] wave;

    showInitial(inr, std::cout);
    */
}
template void Galfit<float>::slit_init(Cube<float> *);
template void Galfit<double>::slit_init(Cube<double> *);


template <class T>
double Galfit<T>::slitfunc(Rings<T> *dring, T *array, int *bhi, int *blo) {


    double minfunc = 0;
    double pixScale = in->Head().PixScale()*arcconv;

    float slitwidth = in->pars().getSlitWidth()/pixScale;
    int bsize[2] = {bhi[0]-blo[0], bhi[1]-blo[1]};
    T *slit = new T[bsize[1]*in->DimZ()];
    int xcenter = lround(dring->xpos.back()-blo[0]);
    int startx = floor(xcenter-slitwidth/2.);
    startx = startx>=0 ? startx : 0.;
    int stopx  = ceil(xcenter+slitwidth/2.);
    stopx = stopx<in->DimX() ? stopx : in->DimX()-1;

    for (int z=0; z<in->DimZ(); z++) {
        for (int y=0; y<bsize[1]; y++) {
            float sum =  0;
            for (int x=startx; x<=stopx; x++) sum += array[x+y*bsize[0]+z*bsize[0]*bsize[1]];
            slit[z+y*in->DimZ()] = sum;
        }
    }

    int numBlanks=0, numPix_tot=0;

    double r1 = dring->radii.front()/pixScale;
    double r2 = dring->radii.back()/pixScale;
    double y0 = dring->ypos.back()-blo[1];

    for (int y=0; y<bsize[0]; y++) {
        bool isIn = (y>=(y0+r1) && y<=(y0+r2)) || (y<=(y0-r1) && y>=(y0-r2));
        if (!isIn) continue;

        float obsSum = 0;
        float modSum = 0;
        float factor = 0;

        for (int zx=0;zx<in->DimZ();zx++) {
            obsSum += (*line_im)(zx,y+blo[1],0);
            modSum += slit[zx+y*in->DimZ()];
        }
        if (modSum!=0) factor = obsSum/modSum;

        for (int zx=0;zx<in->DimZ();zx++) {
            slit[zx+y*in->DimZ()] *= factor;
            T obs = (*line_im)(zx,y+blo[1],0)!=0 ? (*line_im)(zx,y+blo[1],0) : line_im->stat().getSpread();
            T mod = slit[zx+y*in->DimZ()];

            if (obs==0) {
                if (mod==0) continue;
                else numBlanks++;
            }

            numPix_tot++;
            switch(par.FTYPE) {
                case 1:
                    minfunc += std::pow(mod-obs,2)/std::sqrt(obs);
                    break;
                case 2:
                    minfunc += fabs(mod-obs);
                    break;
                case 3:
                    minfunc += fabs(mod-obs)/(mod+obs);
                    break;
                case 4:
                    minfunc += std::pow(mod-obs,2);
                    break;
                default:
                    break;
            }

        }
    }


    delete [] slit;
    return std::pow((1+numBlanks/T(numPix_tot)),par.BWEIGHT)*minfunc/((numPix_tot-numBlanks));
}
template double Galfit<float>::slitfunc(Rings<float>*,float*,int*,int*);
template double Galfit<double>::slitfunc(Rings<double>*,double*,int*,int*);


template <class T>
void Galfit<T>::writeModel_slit() {


    std::string outfold = in->pars().getOutfolder();
    std::string object = in->Head().Name();
    double pixScale = in->Head().PixScale()*arcconv;

    int bhi[2] = {in->DimX(), in->DimY()};
    int blo[2] = {0,0};
    Model::Galmod<T> *mod = getModel(outr,bhi,blo,nullptr,true);

    Cube<T> *modc = mod->Out();
    modc->fitswrite_3d((outfold+object+"_3dmod.fits").c_str());

    float slitwidth = in->pars().getSlitWidth()/pixScale;
    int axis[2] = {in->DimZ(),in->DimY()};
    Image2D<T> *slit = new Image2D<T>(axis);
    line_im->Head().setCrpix(1,findMean(&outr->ypos[0],outr->nr)+1);
    slit->copyHeader(line_im->Head());

    int xcenter = lround(findMean(&outr->xpos[0],outr->nr));
    int startx = floor(xcenter-slitwidth/2.);
    startx = startx>=0 ? startx : 0.;
    int stopx  = ceil(xcenter+slitwidth/2);
    stopx = stopx<in->DimX() ? stopx : in->DimX()-1;

    for (int zx=0; zx<in->DimZ(); zx++) {
        for (int y=0; y<in->DimY(); y++) {
            float sum =  0;
            for (int x=startx; x<=stopx; x++) sum += (*modc)(x,y,zx);
            (*slit)(zx,y) = sum;
        }
    }

    for (int y=0; y<in->DimY(); y++) {
        float obsSum = 0;
        float modSum = 0;
        float factor = 0;
        for (int zx=0;zx<in->DimZ();zx++) {
            obsSum += (*line_im)(zx,y,0);
            modSum += (*slit)(zx,y);
        }
        if (modSum!=0) factor = obsSum/modSum;

        for (int zx=0;zx<in->DimZ();zx++) (*slit)(zx,y) *= factor;

    }
    slit->fitswrite_2d((outfold+object+"mod.fits").c_str());


    std::ofstream outpv_m((outfold+"pv_mod.txt").c_str());
    std::ofstream outpv_o((outfold+"pv.txt").c_str());
    float xmin=1.E10,xmax=0,xmmin=1.E10,xmmax=0;
    for (int y=0; y<slit->DimY() ; y++) {
        for (int x=0;x<slit->DimX();x++) {
            int i = x+y*slit->DimX();
            float xphys = line_im->getXphys(x);
            float yphys = line_im->getYphys(y)*arcconv;
            if (fabs(yphys)<outr->radii[outr->nr-1]+5*outr->radsep ) {
                outpv_m << yphys << "   " << xphys << "  " << (*slit)(i) << endl;
                outpv_o << yphys << "   " << xphys << "  " << (*line_im)(i) << endl;
                if (yphys>xmax) xmax=yphys;
                if (yphys<xmin) xmin=yphys;
            }
        }
        outpv_m << endl;
        outpv_o << endl;
    }
    outpv_m.close();
    outpv_o.close();

    std::ofstream outpv((outfold+"rcpv.txt").c_str());
    for (int i=0; i<outr->nr; i++) {
        float vel1 = (outr->vrot[i]*sin(outr->inc[i]*M_PI/180.))+outr->vsys[i];
        float vel2 = outr->vsys[i]-(outr->vrot[i]*sin(outr->inc[i]*M_PI/180.));
        if (outr->phi.back()<90 || outr->phi.back()>270) std::swap(vel1,vel2);
        float radius = outr->radii[i];
        outpv << -radius << "   " << vel1 << endl;
        outpv <<  radius << "   " << vel2 << endl;
    }
    outpv.close();

    int free[nfree];
    int k;
    for (int nm=0, k=0; nm<9; nm++) {
        if (mpar[nm]) free[k++]=nm;
    }

    const int err_col=13;
    std::ofstream gnu;
    std::string mfile = outfold+"gnuscript.gnu";
    gnu.open(mfile.c_str());

    float xtics = lround(outr->nr/5.)+1;
    xtics *= outr->radsep;
    while (outr->radii.back()/xtics>5) xtics*=2;
    while (outr->radii.back()/xtics<2) xtics/=2;

    /// Setting global option
    gnu << "set terminal postscript eps enhanced color font 'Helvetica,14'" << endl
        << "set output '" << outfold << object << "_vels.eps'" << endl
        << "unset key" << endl
        << "set size 0.60, 1" << endl
        << "set style line 1 lc 9 lt 4 pt 7 lw 1" << endl
        << "set style line 2 lc 8 lt 9 pt 9 lw 1" << endl
        << "set macros" << endl
        << "XTICS   = 'set xtics " << to_string(xtics) << "; set mxtics 2; set format x \"%g\" '" << endl
        << "NOXTICS = 'unset xlabel; set xtics  " << to_string(xtics) << "; set mxtics 2; set format x '' '" << endl
        << "LABELF  = 'set xlabel font \"Helvetica,13\"; "
        <<            "set ylabel font \"Helvetica,13\" '" << endl
        << "TICSF   = 'set xtics font \"Helvetica,12\"; "
        <<            "set ytics font \"Helvetica,12\" '" << endl
        << "TMARGIN = 'set tmargin at screen 0.95; set bmargin at screen 0.47; "
        <<            "set lmargin at screen 0.10; set rmargin at screen 0.50'" << endl
        << "MMARGIN = 'set tmargin at screen 0.47; set bmargin at screen 0.27; "
        <<            "set lmargin at screen 0.10; set rmargin at screen 0.50'" << endl
        << "BMARGIN = 'set tmargin at screen 0.27; set bmargin at screen 0.10; "
        <<            "set lmargin at screen 0.10; set rmargin at screen 0.50'" << endl
        << "set multiplot layout 3,1 rowsfirst" << endl;

    gnu << "@LABELF" << endl << "@TICSF" << endl;

    /// Plotting rotational velocity
    float maxvel = *max_element(&outr->vrot[0], &outr->vrot[0]+outr->nr);
    maxvel += 0.1*maxvel;
    gnu << "@TMARGIN" << endl << "@NOXTICS" << endl
        << "set yrange [-5:" << maxvel << "]" << endl
        << "set ylabel 'V_c  [km/s]'" << endl
        << "set ytics 50" << endl << "set mytics 5" << endl
        << "plot '" << outfold << "rings_final1.txt' ";

    if (par.flagERRORS && mpar[0]) {
        int nc=err_col;
        for (int i=0; i<nfree; i++) if (free[i]==0) nc+=2*i;
        gnu << "u 2:3:($3+$"+to_string(nc)+"):($3+$"+to_string(nc+1)+") w errorbars ls 1, '"
            << outfold <<"rings_final1.txt' u 2:3 w lp ls 1";
    }
    else gnu << "u 2:3 w lp ls 1";

    if (second) {
        gnu << ", '" << outfold << "rings_final2.txt' ";
        if (par.flagERRORS && mpar[0]) {
        gnu << "u 2:3:($3+$13):($3+$14) w errorbars ls 2, '"
            << outfold <<"rings_final2.txt' u 2:3 w lp ls 2";
        }
        else gnu << "u 2:3 w lp ls 2";
    }

    gnu << endl << "set title ''" << endl;
    // Plotting dispersion velocity
    float maxa = *max_element(&outr->vdisp[0], &outr->vdisp[0]+outr->nr);
    maxa += 0.1*maxa;
    gnu << "@MMARGIN" << endl << "@NOXTICS" << endl
        << "set yrange [0:"<<maxa<<"]\n"
        << "set ylabel '{/Symbol s} [km/s]'\n"
        << "set ytics 5" << endl << "set mytics 5" << endl
        << "plot '"<<in->pars().getOutfolder()<<"rings_final1.txt' ";

    if (par.flagERRORS && mpar[1]) {
        int nc=err_col;
        for (int i=0; i<nfree; i++) if (free[i]==1) nc+=2*i;
        gnu << "u 2:4:($4+$"+to_string(nc)+"):($4+$"+to_string(nc+1)+") w errorbars ls 1, '"
            << outfold <<"rings_final1.txt' u 2:4 w lp ls 1";
    }
    else gnu << "u 2:4 w lp ls 1";

    if (second) {
        gnu << ", '" << outfold << "rings_final2.txt' ";
        if (par.flagERRORS && mpar[1]) {
            gnu << "u 2:4:($3+$15):($3+$16) w errorbars ls 2, '"
                << outfold <<"rings_final2.txt' u 2:4 w lp ls 2";
        }
        else gnu << "u 2:4 w lp ls 2";
    }
    gnu << endl;


    // Plotting systemic velocity
    maxa = *max_element(&outr->vsys[0], &outr->vsys[0]+outr->nr);
    maxa += (0.1*maxa+10);
    float mina = *min_element(&outr->vsys[0], &outr->vsys[0]+outr->nr);
    mina -= (0.1*mina+10);
    gnu << "@BMARGIN" << endl << "@XTICS" << endl
        << "set yrange [" << mina << ":" << maxa << "]" << endl
        << "set ylabel 'V_{sys} [km/s]'" << endl
        << "set xlabel 'Radius [arcsec]" << endl
        << "plot '" << outfold << "rings_final1.txt' u 2:12 w lp ls 1";

    if (second)
        gnu << ", '" << outfold << "rings_final2.txt' u 2:12 w lp ls 2";
    gnu << endl;
    gnu << "unset multiplot; reset" << endl;
    gnu.close();

    /// Plotting pv contours
    std::string conlevels;
    float sig;
    if (line_im->pars().getParSE().flagUserGrowthT) sig=line_im->pars().getParSE().threshold;
    else sig = 2.0*line_im->stat().getSpread();
    k=0;
    while (sig<line_im->stat().getMax()) {
        conlevels += to_string(sig)+",";
        sig *= 3;
        k++;
        if (k>10000) break;
    }
    conlevels.erase(conlevels.end()-1);

    float vmin = mina-maxvel-50;
    float vmax = maxa+maxvel+50;
    if (vmin>vmax) std::swap(vmin,vmax);
    mfile=outfold+"pv.gnu";
    outpv.open(mfile.c_str());

    outpv << "set contour base"  << endl
          << "set cntrparam levels discrete "<< conlevels << endl
          << "unset surface"  << endl
          << "set table 'cont.tab'" << endl
          << "splot '"<<outfold+"pv.txt'" << endl
          << "set table 'cont_mod.tab'" << endl
          << "splot '"<<outfold+"pv_mod.txt'" << endl
          << "unset table" << endl
          << "unset key" << endl;

    mfile = outfold+"pv_a_cfr.eps";
    outpv << "set terminal postscript eps enhanced color font 'Helvetica,14'" << endl
          << "set output '" << mfile << "'" << endl
          << "set style line 1 lc rgb '#7F7F7F' lt 2 pt 0 lw 1" << endl
          << "set style line 2 lc rgb '#B22222' lt -1 pt 0 lw 1" << endl
          << "set style line 3 lc rgb '#00008B' lt -1 pt 5 ps 0.6 lw 1" << endl
          << "set xlabel 'Offset [arcsec]'" << endl
          << "set ylabel 'Velocity [km/s]'" << endl
          << "set yrange ["<< to_string(vmin)+":"<< to_string(vmax) << "]" << endl
          << "set yzeroaxis lt -1" << endl
          << "set xrange ["<< to_string(xmin)+":"<< to_string(xmax) << "]" << endl
          << "plot " << to_string(findMean(&outr->vsys[0], outr->nr)) << " ls -1 lc -1, 'cont.tab' w l ls 1, 'cont_mod.tab' w l ls 2, '"<<outfold<<"rcpv.txt' w p ls 3 "<< endl;

#ifdef HAVE_GNUPLOT
    Gnuplot gp;
    gp.begin();
    mfile = "load '"+outfold+"gnuscript.gnu'";
    if (!in->pars().getflagGalMod()) gp.commandln(mfile.c_str());
    mfile = "load '"+outfold+"pv.gnu'";
    gp.commandln(mfile.c_str());
    gp.end();
    remove ("cont.tab");
    remove ("contm.tab");
    remove ("cont_mod.tab");
    remove ("contm_mod.tab");
#endif

    remove ((outfold+"pv.txt").c_str());
    remove ((outfold+"pv_mod.txt").c_str());
    remove ((outfold+"rcpv.txt").c_str());
    remove ((outfold+"pv.gnu").c_str());

    delete slit;

}
template void Galfit<float>::writeModel_slit();
template void Galfit<double>::writeModel_slit();


}

