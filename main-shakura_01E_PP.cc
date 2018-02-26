/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 * Post-processing main function for analyzing thin disk simulations in
 * Cylindrical KS metric
 *++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include "../exe-dir/header.hh"
#include "../PostProcess/PostProcess.hh"
#include "../PostProcess/MetricPP.hh"

using namespace std;

void comovingTensor(Vector3d, double, Tensor4d, Vector4d,
                    vector<vector<double> >, Tensor3d &);

void  setZoneConversion(vector<int> &mzone,
                        int  nx, int  ny, int  nz,
                        int  px, int  py, int  pz,
                        int opx, int opy, int opz);

int  main(int argc, char *argv[])
{
    /* global grid parameters
     *   nx, ny, & nz are the number of zones in each dimension*/
    int nb = 1;
    int nx = 256;
    int ny = 192;
    int nz = 1;
    int px = atoi( argv[1] );
    int py = atoi( argv[2] );
    int pz = atoi( argv[3] );
    int nblock = atoi( argv[4] );
    
    /* create grid to match original run*/
    int opx = 1;
    int opy = 1;
    int opz = 1;
    int nlevels     = 1;
    int nprerefine  = 0;
    
    int ndimensions = 1;
    if(ny > 1) ndimensions = 2;
    if(nz > 1) ndimensions = 3;
    
    int iRestart      = 0;
    int nStart        = 0;
    int iStepSize     = 1;
    int nStop         = 5000;
    int nStartAverage = 0; /* indicates at which step averaging begins*/
    int iFluxShell    = 200;
    
    bool isVerboseOn         = false;
    bool isMagFieldOn        = false;
    bool isMagFieldPrimitive = true;
    bool isGasIdeal          = true;
    bool isPressureDumped    = true;
    bool isTemperatureDumped = true;
    bool isRadOn             = true;
    bool isComptonOn         = true;
    
    double pi        = 2.0*asin(1.0);
    double msolar    = 1.989e+33;
    double bhmass    = 6.62*msolar;
    double Ledd      = 1.3e+38*bhmass/msolar;
    double tol       = 1.e-5;
    double tiny      = Constants::tinyNumber;
    double mu        = 1.69;
    double mH        = Constants::protonMass;
    double boltz     = Constants::boltzmanConstant;
    double gee       = Constants::gravityConstant;
    double aR        = Constants::radiationConstant;
    double kb        = Constants::boltzmanConstant;
    double aCon      = Constants::radiationConstant;
    double cee       = Constants::speedOfLight;
    double m_e       = Constants::electronMass;
    double lunit     = gee*bhmass/(cee*cee);
    double tunit     = lunit/cee;
    double dunit     = (cee*cee)/(lunit*lunit*gee);
    double munit     = dunit*lunit*lunit*lunit;
    double edunit    = dunit*(lunit*lunit)/(tunit*tunit);
    double tempunit  = mH*cee*cee/boltz;
    m_e  /= munit;
    cee  *= (tunit/lunit);
    kb   *= (tempunit*tunit*tunit)/(munit*lunit*lunit);
    aCon *= (lunit*tunit*tunit*pow(tempunit,4))/munit;
    
    Constants::timeUnit     = tunit;
    Constants::massUnit     = munit;
    Constants::lengthUnit   = lunit;
    Constants::temperatureUnit = tempunit;
    
    bhmass /= munit;
    double gamma     = 5.0/3.0;
    double bhspin    = 0.0;
    
    /*Thin disk parameters*/
    double rho0      = 1.0e-4;
    double h2        = 0.4;
    double h1        = h2;
    double hthick    = 20.*h2;
    double rintd     = 6.0;
    double rmintd    = 5.0;
    double rmaxtd    = 20.0;
    double alpha     = 0.02;
    
    // radiation parameters (cgs units)
    // Scattering opacity
    double Ks       = 0.34; 		// cm^2/g
    // Planck opacity
    double PlCoef   = 6.4e+22;       	// cm^2/g
    double PlTexp   = -7./2.; 		// bremsstrahlung = -7/2
    double PlDexp   = 1.; 		// bremsstrahlung = 1
    // Rosseland mean opacity
    double RoCoef   = 1.6e+21;       	// cm^2/g
    double RoTexp   = -7./2.; 		// bremsstrahlung = -7/2
    double RoDexp   = 1.; 		// bremsstrahlung = 1
    
    //grid parameters
    double eta0      = 1.0;
    double etaz0     = 0.0;
    double r1        = 0.0;
    double rI        = 1.0;
    int ilogr        = 1;
    int ilogz        = 1;
    int    AngTrans  = 0;
    double AngR0     = 1.;
    double AngR1     = 1.;
    double AngR2     = 0.;
    double AngPar    = 0.1;
    double AngCut    = 0.;
    double AngIndex  = 0.;
    
    double radCut    = 20.;
    double radAccel  = 10.;
    
    /* Initialization functions for grid*/
    MPIClass mpi;
    mpi.startMPI(argc, argv, px, py, pz);
    int rank   = mpi.getProcID();
    int nprocs = mpi.getNumProcs();
    
    /*grid options*/
    double ratio       = 1.0;
    
    /*determine black hole radius*/
    double bhradius = bhmass + sqrt(bhmass*bhmass - bhspin*bhspin + tiny);
    if(rank == 0) cout << "R_BH     = " << bhradius << endl;
    
    double rmin    = rmintd;
    double rmax    = rmaxtd;
    double r0      = bhradius;
    
    /*simulation box parameters*/
    
    double boxx = rmax - rmin;
    double boxy = 0.3*pi;
    double boxz = 0.25*pi;
    if(nz == 1) boxz = 1.;
    
    double xmin = rmin;
    double ymin = 0.5*(pi - boxy);
    double zmin = -0.5*boxz;
    
    
    if(ilogr == 1) {
        boxx = log(rmax/r0) - log(rmin/r0);
        xmin = eta0 + log(xmin/r0);
    }
    
    
    double xmax = xmin + boxx;
    double ymax = ymin + boxy;
    double zmax = zmin + boxz;
    
    if(rank == 0) cout << "Setting up the mesh..." << endl;
    KLMMesh *mesh;
    if(ndimensions == 1) {
        mesh = new KLMMesh(mpi,nlevels,nb,nx,xmin,xmax,0);
    } else if(ndimensions == 2) {
        mesh = new KLMMesh(mpi,nlevels,nb,nx,ny,xmin,xmax,ymin,ymax,0);
    } else {
        mesh = new KLMMesh(mpi,nlevels,nb,nx,ny,nz,
                           xmin,xmax,ymin,ymax,zmin,zmax,0);
    }
    
    Field field(*mesh, mpi);
    
    EOS *eos;
    if(isGasIdeal) {
        eos = new IdealGas(*mesh, field, mpi);
    } else {
        eos = new RelativisticGas(*mesh, field, mpi);
    }
    (*eos).setGamma(gamma);
    
    // opacities are set in cgs units (cm^2/gm)
    OpacPow opacPlanck;
    opacPlanck.setOpacityConstant(PlCoef);
    opacPlanck.setTemperatureExponent(PlTexp);
    opacPlanck.setDensityExponent(PlDexp);
    
    OpacPow opacRosseland;
    opacRosseland.setOpacityConstant(RoCoef);
    opacRosseland.setTemperatureExponent(RoTexp);
    opacRosseland.setDensityExponent(RoDexp);
    
    OpacCon opaccon;
    opaccon.setOpacity(Ks);
    
    Opacity *opacity;
    opacity = new Opacity();
    opacity->addAbsorption("planckMean", &opacPlanck);
    opacity->addAbsorption("rossMean", &opacRosseland);
    opacity->addScattering(&opaccon);
    
    Metric *metric;
    SphericalKS metricks (*mesh, field, mpi);
    SphericalBL metricbl(*mesh, field, mpi);
    metricbl.setBHSpin(bhspin);
    metricbl.setBHMass(bhmass);
    metricks.setBHSpin(bhspin);
    metricks.setBHMass(bhmass);
    metricks.setR0(r0);
    metricks.setEta0(eta0);
    metricks.setR1(r1);
    metricks.setRindex(rI);
    if(ilogr == 1) metricks.turnOnLogRCoordinate();
    
    metricks.setAngularTransform(AngTrans);
    metricks.setAngR0(AngR0);
    metricks.setAngR1(AngR1);
    metricks.setAngR2(AngR2);
    metricks.setAngularParameter(AngPar);
    metricks.setAngularCut(AngCut);
    metricks.setAngularIndex(AngIndex);
    
    metric = &metricks;
    metric->setMetric();
    
    vector<vector<vector<int> > > oppZoneIndex = field.getOppZoneIndex();
    /*Define post-processor object, read Master file, declare fields.
     *This part should only be done once.*/
    string baseDirName = string("/Users/mbhupe/JILA_projects/visc_test/");
    string fileStem = "Shakura_01E";
    /*create post-processing object*/
    
    PPCLASS   pp(mpi);
    pp.setVerbose(isVerboseOn);
    pp.readMasterFile(baseDirName,fileStem);
    vector<double> rho, energy,radiationEnergy, pressure, temperature, magpress, boost, sqrtg,
    velocity, velocity3, magfield, magfield3, radvelocity, radvelocity3;
    int totalSteps = pp.getNumTimeSteps();
    if(totalSteps < nStop) {
        if(mpi.getProcID() == 0) {
            cout << "WARNING: nStop exceeds value of totalSteps." << endl;
            cout << "Defaulting nStop to: " << totalSteps << endl;
        }
        nStop = totalSteps;
    }
    
    double timeold = -1.0e+20;
    double timemax = pp.getCurrentTime(totalSteps-1);
    
    int nrad   = nx;
    int ntheta = ny;
    int nphi   = nz;
    if(nprerefine > 0) {
        nrad   *= 2*nprerefine;
        ntheta *= 2*nprerefine;
        nphi   *= 2*nprerefine;
    }
    vector<double> rad(nrad),theta(ntheta),photoZmin(nrad),photoZmax(nrad);
    
    pp.readData(nStart, px, py, pz, opx, opy, opz);
    pp.getScalarField("massDensity", rho);
    vector<int> mzone(rho.size());
    if(nlevels == 1) {
        setZoneConversion(mzone, nx, ny, nz, px, py, pz, opx, opy, opz);
    } else {
        for(int n = 0; n < rho.size(); n++) {
            mzone[n] = n;
        }
    }
    int nFaces = (*mesh).getZone(0)->getNumberOfFaces();
    vector<Vector3d> faceArea(nFaces);
    
    vector<Vector3d> posXYZ(rho.size()), posRPZ(rho.size()), posRTP(rho.size());
    for(int izone = 0; izone < rho.size(); izone++) {
        Vector3d pos3d;
        pos3d = (*mesh).getZone(mzone[izone])->getPosition();
        (*metric).getCartesianPosition(pos3d);
        posXYZ[izone] = pos3d;
        
        pos3d = (*mesh).getZone(mzone[izone])->getPosition();
        (*metric).getPhysicalPosition(pos3d);
        posRPZ[izone] = pos3d;
        
        pos3d = (*mesh).getZone(mzone[izone])->getPosition();
        (*metric).getSphericalPosition(pos3d);
        posRTP[izone] = pos3d;
    }
    
    double dr, dr0, x0, dtheta, y0;
    if(ratio == 1.0) {
        dr0 = (xmax - xmin)/nrad;
    } else {
        dr0 = (xmax-xmin)*(1.0-ratio)/(1.0-pow(ratio,nrad));
    }
    x0 = xmin;
    dtheta = (ymax - ymin)/ntheta;
    for(int irad = 0; irad < nrad; irad++) {
        dr  = dr0*pow(ratio, irad);
        rad[irad] = x0 + 0.5*dr;
        if(ilogr) rad[irad] = r1 + r0*exp(pow((rad[irad]-eta0),rI));
        y0 = ymin;
        for(int j = 0; j < ntheta; j++) {
            double a11 = 0.5+1./pi*atan((rad[irad]-AngR2)/AngR1);
            double a12 = 2.-(2.-AngPar)*pow(rad[irad]/AngR0,-AngIndex*a11);
            double x2  = y0 + 0.5*dtheta;
            theta[j] = x2 + 0.5*(1.0-a12)*sin(2.0*x2);
            y0 += dtheta;
        }
        x0 += dr;
    }
    
    /*Time & shell-averaged quantities*/
    int nStopOld, nSkipOld, nAverageOld;
    int nSkip  = 0;
    int nAverage = 0;
    vector<double> tsAvgQ_MRI2(nrad, 0.0);
    vector<double> tsAvgQ_MRI3(nrad, 0.0);
    vector<double> tsAvgCoolingRate(nrad, 0.0);
    vector<double> tsAvgalpha1(nrad, 0.0);
    vector<double> tAvgMassFlux(nrad, 0.0);
    vector<double> tAvgPhotoZmin(nrad, 0.0);
    vector<double> tAvgPhotoZmax(nrad, 0.0);
    vector<double> tAvgPhotoZmin_Scat(nrad, 0.0);
    vector<double> tAvgPhotoZmax_Scat(nrad, 0.0);
    vector<double> tAvgPhotoZmin_Abs(nrad, 0.0);
    vector<double> tAvgPhotoZmax_Abs(nrad, 0.0);
    vector<double> tAvgPhotoZmax_Eff(nrad, 0.0);

    vector<double> Ptot0(nrad,0.0);
    vector<double> Pmag0(nrad,0.0);
    vector<double> radV0(nrad,0.0);
    vector<double> vz0(nrad,0.0);
    vector<double> tzAvgRho(ntheta, 0.0);
    vector<double> tzAvgarad(ntheta, 0.0);
    vector<double> tzAvgagas(ntheta, 0.0);
    vector<double> tzAvgamag(ntheta, 0.0);
    vector<double> faceArea4(nrad, 0.0);
    vector<vector<double> > tsAvgRho(nrad);
    vector<vector<double> > tsAvgErad(nrad);
    vector<vector<double> > tsAvgPratio(nrad);
    vector<vector<double> > tsAvgv1(nrad);
    vector<vector<double> > tsAvgv3(nrad);
    vector<vector<double> > tsAvgb1(nrad);
    vector<vector<double> > tsAvgb3(nrad);
    vector<vector<double> > tsAvgvrad1(nrad);
    vector<vector<double> > tsAvgvrad3(nrad);
    vector<vector<double> > tsAvgRadV(nrad);
    vector<vector<double> > tsAvgVz(nrad);   
 
    for(int irad=0;irad < nrad;irad++){
        tsAvgRho[irad].resize(ntheta,0.0);
        tsAvgErad[irad].resize(ntheta,0.0);
        tsAvgPratio[irad].resize(ntheta,0.0);
        tsAvgv1[irad].resize(ntheta,0.0);
        tsAvgv3[irad].resize(ntheta,0.0);
        tsAvgb1[irad].resize(ntheta,0.0);
        tsAvgb3[irad].resize(ntheta,0.0);
        tsAvgvrad1[irad].resize(ntheta,0.0);
        tsAvgvrad3[irad].resize(ntheta,0.0);
        tsAvgRadV[irad].resize(ntheta,0.0);
        tsAvgVz[irad].resize(ntheta,0.0);
    }
    
    if (nStart > nStop){
        cout << "ERROR:  nStop has been set lower than nStart." << endl;
        exit(0);
    }
    
    for(int n = nStart; n < nStop; n+=iStepSize) {
        double time = pp.getCurrentTime(n);
        if( (time-timeold) < tol*timemax ) {
            nSkip ++;
            if(mpi.getProcID() == 0) cout << "skipping time step = " << n << endl;
        } else {
            if(mpi.getProcID() == 0) cout << "time step = " << n << endl;
            pp.readData(n, px, py, pz, opx, opy, opz);
            pp.getScalarField("massDensity", rho);
            pp.getScalarField("boostFactor", boost);
            pp.getScalarField("squareRootMetricDet", sqrtg);
            pp.getVectorField("velocity", velocity);
            if(nz == 1) pp.getScalarField("velocity3", velocity3);
            if(isPressureDumped) {
                pp.getScalarField("pressure", pressure);
                energy.resize( pressure.size() );
            } else {
                pp.getScalarField("energyDensity", energy);
                pressure.resize( energy.size() );
            }
            if(isTemperatureDumped)
                pp.getScalarField("temperature", temperature);
            else
                temperature.resize( energy.size() );
            if(isRadOn) {
                pp.getScalarField("radiationEnergy",radiationEnergy);
                pp.getVectorField("radiationVelocity", radvelocity);
                if(nz == 1) pp.getScalarField("radiationVelocity3", radvelocity3);
            } else {
                radiationEnergy.resize( rho.size(), 0. );
                if(nz == 1) {
                    radvelocity.resize( 2*rho.size(), 0. );
                    radvelocity3.resize( rho.size(), 0. );
                } else {
                    radvelocity.resize( 3*rho.size(), 0. );
                }
            }
            if(isMagFieldOn) {
                pp.getScalarField("magneticPressure", magpress);
                pp.getVectorField("magneticField", magfield);
                if(nz == 1) pp.getScalarField("magneticField3", magfield3);
            } else {
                magpress.resize( rho.size(), 0. );
                if(nz == 1) {
                    magfield.resize( 2*rho.size(), 0. );
                    magfield3.resize( rho.size(), 0. );
                } else {
                    magfield.resize( 3*rho.size(), 0. );
                }
            }
            
            // fill in EOS variables
            for(int izone = 0; izone < rho.size(); izone++) {
                if(!isGasIdeal) {
                    double eps = energy[izone]/rho[izone];
                    gamma = (5.+4.*eps)/3./(1.+eps);
                }
                if(isPressureDumped) {
                    energy[izone] = pressure[izone]/(gamma-1.);
                } else {
                    pressure[izone] = energy[izone]*(gamma-1.);
                }
                if(!isTemperatureDumped) {
                    double nt = rho[izone]*dunit/(mH*mu);
                    temperature[izone] = (gamma-1.)*energy[izone]*edunit/(boltz*nt)/tempunit;
                }
            }
            
            //convert vectors to spherical coords
            vector<Vector3d> velocityRTP(rho.size());
            vector<Vector3d> velocityXYZ(rho.size());
            vector<Vector3d> magfieldRTP(rho.size());
            vector<Vector3d> magfieldXYZ(rho.size());
            vector<Vector3d> radVelRTP(rho.size());
            vector<Vector3d> radVelXYZ(rho.size());
            
            for(int izone = 0; izone < rho.size(); izone++) {
                double x  = posXYZ[izone].getX();
                double y  = posXYZ[izone].getY();
                double z  = posXYZ[izone].getZ();
                double vx,vy,vz,vrx,vry,vrz;
                
                if(ndimensions == 2) {
                    vx = velocity[2*izone+0];
                    vy = velocity3[izone];
                    vz = velocity[2*izone+1];
                    vrx = radvelocity[2*izone+0];
                    vry = radvelocity3[izone];
                    vrz = radvelocity[2*izone+1];
                } else {
                    vx = velocity[3*izone+0];
                    vy = velocity[3*izone+1];
                    vz = velocity[3*izone+2];
                    vrx = radvelocity[3*izone+0];
                    vry = radvelocity[3*izone+1];
                    vrz = radvelocity[3*izone+2];
                }
                
                double r = posRTP[izone].getX();
                double sintheta = sin(posRTP[izone].getY());
                double costheta = cos(posRTP[izone].getY());
                double sinphi = sin(posRTP[izone].getZ());
                double cosphi = cos(posRTP[izone].getZ());
                double rcyl2 = x*x + y*y;
                double drdx = x/r;
                double drdy = y/r;
                double drdz = z/r;
                double dtdx = (x/(r*r)-rcyl2/(r*r*r)*drdx)/(sintheta*costheta);
                double dtdy = (y/(r*r)-rcyl2/(r*r*r)*drdy)/(sintheta*costheta);
                double dtdz = -(1.0-z/r*drdz)/(r*sintheta);
                double dpdx;
                if(ndimensions == 2) {
                    dpdx = 0.0;
                } else {
                    dpdx = -(1.0-x*x/rcyl2)/(sqrt(rcyl2)*sinphi);
                }
                double dpdy =  (1.0-y*y/rcyl2)/(sqrt(rcyl2)*cosphi);
                double dpdz = 0.0;
                // spherical-polar components in grid frame
                double vr  = drdx*vx+drdy*vy+drdz*vz;
                double vth = dtdx*vx+dtdy*vy+dtdz*vz;
                double vph = dpdx*vx+dpdy*vy+dpdz*vz;
                
                double vrr  = drdx*vrx+drdy*vry+drdz*vrz;
                double vrth = dtdx*vrx+dtdy*vry+dtdz*vrz;
                double vrph = dpdx*vrx+dpdy*vry+dpdz*vrz;
                velocityRTP[izone].setX(vr);
                velocityRTP[izone].setY(vth);
                velocityRTP[izone].setZ(vph);
                
                radVelRTP[izone].setX(vrr);
                radVelRTP[izone].setY(vrth);
                radVelRTP[izone].setZ(vrph);
            }
            
            //reset metric components to SphericalKS
            SphericalKS metricks2(*mesh, field, mpi);
            metricks2.setBHSpin(bhspin);
            metricks2.setBHMass(bhmass);
            //metricks2.setBHTilt(bhtilt);
            
            if(n >= nStartAverage) nAverage += 1;
            
            //calculate gradients for acceleration terms
            int numZones = (*mesh).getNumZones();
            vector<double> Pgas(numZones, 0.), Pmag(numZones, 0.);
            for(int izone = 0; izone < rho.size(); izone++) {
                Pgas[mzone[izone]] = pressure[izone];
                Pmag[mzone[izone]] = magpress[izone];
            }
            vector<Vector3d> gradPgas,gradPmag;
            gradPgas = field.getVectorGrad(Pgas);
            gradPmag = field.getVectorGrad(Pmag);
            
            vector<vector<double> > kappaTRho(nrad);
            vector<vector<double> > kappaARho(nrad);
            vector<vector<double> > kappaScat(nrad);
            vector<vector<double> > shellHeating(nrad);
            vector<vector<double> > shellStress(nrad);
            vector<vector<double> > shellFluxR(nrad);
            vector<vector<double> > shellFluxZ(nrad);
            vector<vector<double> > RadV(nrad);
            vector<vector<double> > Vz(nrad);


            for(int irad=0;irad < nrad;irad++){
                kappaTRho[irad].resize(ntheta,0.0);
                kappaARho[irad].resize(ntheta,0.0);
                kappaScat[irad].resize(ntheta,0.0);
                shellHeating[irad].resize(ntheta,0.0);
                shellStress[irad].resize(ntheta,0.0);
                shellFluxR[irad].resize(ntheta,0.0);
                shellFluxZ[irad].resize(ntheta,0.0);
                RadV[irad].resize(ntheta,0.0);
                Vz[irad].resize(ntheta,0.0);
            }
            
            
            /**********************************************************/
            /* Start shell-averaging */
            /**********************************************************/
            double mass = 0.;
            double luminosity = 0.;
            
            for(int irad = 0; irad < nrad; irad++) {
                double shellAvgScaleHeight2= 0.0;
                double shellAvgRho2        = 0.0;
                double shellAvgRho         = 0.0;
                double shellAvgPress       = 0.0;
                double shellAvgPrad        = 0.0;
                double shellAvgPrad_Pgas   = 0.0;
                double shellAvgTgas        = 0.0;
                double shellAvgTgasZ0      = 0.0;
                double shellAvgTrad        = 0.0;
                double shellAvgTradZ0      = 0.0;
                double shellAvgVphi        = 0.0;
                double shellAvgQ_MRI2      = 0.0;
                double shellAvgQ_MRI3      = 0.0;
                double shellAvgAlfvenSpeed = 0.0;
                double shellAvgAlfvenSpeed2= 0.0;
                double shellAvgAlfvenSpeed3= 0.0;
                double shellAvgalpha1      = 0.0;
                double shellMassFlux       = 0.0;
		        double shellAvgvzz         = 0.0;
		        double shellAvgradV        = 0.0;
                double shellAvgErad        = 0.0;
                double shellAvgPtot0       = 0.0;
		        double shellAvgRho0        = 0.0;
                double shellAvgPmag0       = 0.0;

                vector<double> azAvgRho(ntheta, 0.0);
                vector<double> azAvgErad(ntheta, 0.0);
                vector<double> azAvgPratio(ntheta, 0.0);
                vector<double> azAvgv1(ntheta, 0.0);
                vector<double> azAvgv3(ntheta, 0.0);
                vector<double> azAvgb1(ntheta, 0.0);
                vector<double> azAvgb3(ntheta, 0.0);
                vector<double> azAvgvrad1(ntheta, 0.0);
                vector<double> azAvgvrad3(ntheta, 0.0);
                
                for(int izone = 0; izone < rho.size(); izone++) {
                    if(fabs(posRTP[izone].getX() - rad[irad])/rad[irad] < tol) {
                        mass += boost[izone]*rho[izone]*(*mesh).getZone(mzone[izone])->getVolume();
                        double r    = posRTP[izone].getX();
                        double th   = posRTP[izone].getY();
                        double phi  = posRTP[izone].getZ(); 
                        double vol  = (*mesh).getZone(mzone[izone])->getVolume();
                        Tensor4d metks, metbl, invMet;
                        faceArea = (*mesh).getZone(mzone[izone])->getFaceArea();
                        metricks2.getMetric(posRTP[izone], metks);
                        metricbl.getMetric(posRTP[izone], metbl);
                        invMet = metks.getInverse();
                        shellAvgScaleHeight2 += rho[izone]*rho[izone]*sqrtg[izone]*vol*(posRTP[izone].getY() - pi/2.)*(posRTP[izone].getY() - pi/2.);
                        shellAvgRho          += rho[izone]*sqrtg[izone]*vol;
                        shellAvgRho2         += rho[izone]*rho[izone]*sqrtg[izone]*vol;
                        shellAvgPress        += pressure[izone]*rho[izone]*sqrtg[izone]*vol;
                        shellAvgTgas         += temperature[izone]*rho[izone]*rho[izone]*sqrtg[izone]*vol;
                        shellAvgVphi         += velocityRTP[izone].getZ()*rho[izone]*sqrtg[izone]*vol;
                        
                        /*if(posRTP[izone].getY() > -0.4 && posRTP[izone].getY() < 0.4) {*/
                        double mdot           = boost[izone]*rho[izone]*velocityRTP[izone].dot(faceArea[0]);
			// vzz                  += rho[izone]*velocity[2*izone+1]*sqrtg[izone]*vol;
			// radV                 += rho[izone]*velocityRTP[izone].getX()*sqrtg[izone]*vol;
			
                        if(ilogr == 1) { //To adjust Vr to Veta to match units of boost and face area
                            double r = posRTP[izone].getX();
                            double detadr = pow(log((r-r1)/r0),1./rI-1.)/rI/(r-r1);
                            //cout << "detadr = " << detadr << endl;
                            mdot *= detadr;
                        }
                        if(nphi == 1) mdot *= 2.*pi;
                        shellMassFlux        += mdot;
                        
                        //}
                        if(isRadOn) {
                            shellAvgErad       += radiationEnergy[izone]*rho[izone]*sqrtg[izone]*vol;
                            shellAvgPrad       += radiationEnergy[izone]*rho[izone]*sqrtg[izone]*vol/3.;
                            shellAvgPrad_Pgas  += sqrtg[izone]*vol*rho[izone]*(radiationEnergy[izone]/3.)/pressure[izone];
                            shellAvgTrad       += pow(radiationEnergy[izone]/aR,0.25)*rho[izone]*sqrtg[izone]*vol;
                            
                            if(posRTP[izone].getY() <= 0.01 && posRTP[izone].getY() > -0.01  ) {
                                shellAvgTgasZ0         += temperature[izone];
                                shellAvgTradZ0         += pow(radiationEnergy[izone]/aR,0.25);
                            }
                        }
                        
                        double velr = velocityRTP[izone].getX();
                        double velt = velocityRTP[izone].getY();
                        double velp = velocityRTP[izone].getZ();
                        double enthalpy = 1.0 + gamma*energy[izone]/rho[izone];
                        if(isMagFieldOn)
                        {
                            double magr = magfieldRTP[izone].getX();
                            double magt = magfieldRTP[izone].getY();
                            double magp = magfieldRTP[izone].getZ();
                            double numerator = metks(0,1)*magr + metks(0,2)*magt + metks(0,3)*magp
                            + metks(1,1)*magr*velr + metks(2,2)*magt*velt + metks(3,3)*magp*velp
                            + metks(1,2)*(magr*velt + magt*velr)
                            + metks(1,3)*(magr*velp + magp*velr)
                            + metks(2,3)*(magt*velp + magp*velt);
                            double denominator = metks(0,0) + metks(0,1)*velr +
                            metks(0,2)*velt + metks(0,3)*velp;
                            double mag0up = -numerator/denominator;
                            Vector4d b      = Vector4d( mag0up, magr, magt, magp );
                            double enthalpy = 1.0 + gamma*energy[izone]/rho[izone];
                            double va       = sqrt(2.0*magpress[izone]/(rho[izone]*enthalpy + 2.0*magpress[izone]));
                            shellAvgAlfvenSpeed += va*rho[izone]*sqrtg[izone]*vol;
                            Vector4d b_dn   = metks.dot(b);
                            faceArea        = (*mesh).getZone(mzone[izone])->getFaceArea();
                            double deltaT = (*mesh).getZone(mzone[izone])->getVolume()/faceArea[2].magnitude();
                            double dthetadx2;
                            double x2 = (*mesh).getZone(mzone[izone])->getPosition().getY();
                            deltaT *= dthetadx2;
                            double va2 = sqrt(abs(b_dn[2]*b[2])/(rho[izone]*enthalpy +
                                                                 2.0*magpress[izone]));
                            shellAvgAlfvenSpeed2 += va2*rho[izone]*vol;
                            shellAvgQ_MRI2 += 2*pi*va2/abs(velp)/(rad[irad]*deltaT)*rho[izone]*vol;
                            if(ndimensions == 3) {
                                double deltaP = (*mesh).getZone(mzone[izone])->getVolume()/faceArea[4].magnitude();
                                deltaP *= sin(posRTP[izone].getY());
                                double va3 = sqrt(abs(b_dn[3]*b[3])/(rho[izone]*enthalpy +
                                                                     2.0*magpress[izone]));
                                shellAvgAlfvenSpeed3 += va3*rho[izone]*vol;
                                shellAvgQ_MRI3 += 2*pi*va3/abs(velp)/(rad[irad]*deltaP)*rho[izone]*vol;
                            }
                        }
                        
                        /***********************************************/
                        vector<vector<double> > stress(4),covStress(4),gmn(4),gmnup(4);
                        //four velocity in Cylindrical coordinates
                        Vector4d u;
                        u[0] = boost[izone]/sqrtg[izone];
                        u[1] = u[0]*velr;
                        u[2] = u[0]*velt;
                        u[3] = u[0]*velp;
                        
                        for(int i = 0; i < 4; i++) {
                            gmn[i].resize(4, 0.);
                            gmnup[i].resize(4, 0.);
                            stress[i].resize(4, 0.);
                            covStress[i].resize(4, 0.);
                            for(int j = 0; j < 4; j++) {
                                if(j > i) {
                                    gmn[i][j] = metks(i,j);
                                    gmnup[i][j] = invMet(i,j);
                                } else {
                                    gmn[i][j] = metks(j,i);
                                    gmnup[i][j] = invMet(j,i);
                                }
                            }
                        }
                        double totalPressure = pressure[izone] + magpress[izone];
                        // Calculate stress-energy tensor
                        for(int i = 0; i < 4; i++) {
                            for(int j = 0; j < 4; j++) {
                                stress[i][j] = (rho[izone]*enthalpy)*u[i]*u[j] + totalPressure*gmnup[i][j];
                            }
                        }
                        // Lower indices on stress tensor
                        for(int i = 0; i < 4; i++) {
                            for(int j = 0; j < 4; j++) {
                                for(int k = 0; k < 4; k++) {
                                    for(int l = 0; l < 4; l++) {
                                        covStress[i][j] += gmn[i][k]*gmn[j][l]*stress[k][l];
                                    }}}}
                        // Calculate stress tensor in co-moving frame
                        Tensor3d W = Tensor3d(0.);
                        comovingTensor(posRTP[izone], bhspin, metbl, u, covStress, W);
                        
                        if(isRadOn) totalPressure += radiationEnergy[izone]/3.;
                        shellAvgalpha1 += (W(0,2)/totalPressure)*rho[izone]*sqrtg[izone]*vol;
                        
                        /***********************************************/
                        for(int itheta = 0 ; itheta < ntheta ; itheta++) {
                            if(fabs(posRTP[izone].getY() - theta[itheta])/fabs(theta[itheta]) < tol) {
                                azAvgRho[itheta] += rho[izone];
                                azAvgv1[itheta] += velr;
                                azAvgv3[itheta] += velp;
                                if(itheta < ntheta/2) dtheta = theta[itheta+1]-theta[itheta];
                                else dtheta = theta[itheta]-theta[itheta-1];
                                Vector3d dxCov = (*mesh).getZone(mzone[izone])->getVectorCovariantZoneLength();
                                kappaTRho[irad][itheta] += (boost[izone]/sqrtg[izone])*rho[izone]*dunit*(Ks + PlCoef*pow(temperature[izone]*tempunit,PlTexp)*pow(rho[izone]*dunit,PlDexp))*r*dtheta*lunit;
                                kappaScat[irad][itheta] += (boost[izone]/sqrtg[izone])*rho[izone]*dunit*Ks*r*dtheta*lunit;
                                kappaARho[irad][itheta] += (boost[izone]/sqrtg[izone])*rho[izone]*dunit*(PlCoef*pow(temperature[izone]*tempunit,PlTexp)*pow(rho[izone]*dunit,PlDexp))*r*dtheta*lunit;
                                RadV[irad][itheta]      += velocityRTP[izone].getX()*cos(2.*pi*0.002*time);
                                Vz[irad][itheta]        += velocity[2*izone+1]*cos(2.*pi*0.002*time);

                                if(isMagFieldOn) {
                                    shellStress[irad][itheta] += W(0,2)*r*dtheta;
                                    shellHeating[irad][itheta] += 1.5*velp*W(0,2)*r*dtheta;
                                } else {
                                    shellStress[irad][itheta] += totalPressure*alpha*r*dtheta;
                                    shellHeating[irad][itheta] += 1.5*velp*totalPressure*alpha*r*dtheta;
                                }
                                double vrr      = radVelRTP[izone].getX();
                                double vrth     = radVelRTP[izone].getY();
                                double vrph     = radVelRTP[izone].getZ();
                                double ur0      = sqrt(-1./(metks(0,0)+2.*(metks(0,1)*vrr+metks(0,2)*vrth+metks(0,3)*vrph)+2.*(metks(1,2)*vrr*vrth+metks(2,3)*vrth*vrph+metks(1,3)*vrr*vrph)+metks(1,1)*vrr*vrr+metks(2,2)*vrth*vrth+metks(3,3)*vrph*vrph));
                                double ur_0      = metks(0,0) + metks(0,1)*vrr + metks(0,2)*vrth + metks(0,3)*vrph;
                                ur_0 *= ur0;
                                shellFluxR[irad][itheta] -= 4./3.*radiationEnergy[izone]*vrr*ur0*ur_0;
                                shellFluxZ[irad][itheta] -= 4./3.*radiationEnergy[izone]*r*vrth*ur0*ur_0;
                                if(itheta == ntheta/2-1 || itheta == ntheta/2) {
                                    shellAvgPtot0 += 0.5*(pressure[izone]+radiationEnergy[izone]/3.0 + magpress[izone]);
                                    shellAvgPmag0 += 0.5*magpress[izone];
                                    shellAvgvzz   += 0.5*velocity[2*izone+1];
                                    shellAvgradV  += 0.5*velocityRTP[izone].getX();
                                    faceArea4[irad] += 0.5*sqrt(gmn[1][1]*gmn[2][2])*faceArea[4].magnitude();
				                    shellAvgRho0    += rho[izone];
                                }
                                azAvgErad[itheta] += radiationEnergy[izone];
                                azAvgPratio[itheta] += magpress[izone]/(pressure[izone]+radiationEnergy[izone]/3.0);
                                azAvgvrad1[itheta] += vrr;
                                azAvgvrad3[itheta] += vrth;
                                if(isRadOn) {
                                    if(irad < iFluxShell && (itheta == 0 || itheta == ntheta-1) ||
                                       irad == iFluxShell) {
                                        // NOTE: This takes advantage of the fact that the Schwarzschild
                                        // metric is diagonal.  Needs to be redone if spin is included.
                                        double g = metks.getSquareRootDet();
                                        double Rrt = 4./3.*radiationEnergy[izone]*vrr*ur0*ur_0;
                                        double Rtht = 4./3.*radiationEnergy[izone]*vrth*ur0*ur_0;
                                        double dr = (xmax-xmin)/nx;
                                        if(ilogr > 0) dr *= rad[irad];
                                        double dphi = boxz/nz;
                                        if(ndimensions == 2) dphi = 2.*pi;
                                        if(irad == iFluxShell) luminosity -= g*Rrt*dtheta*dphi;
                                        if(itheta == 0) luminosity += g*Rtht*dr*dphi;
                                        if(itheta == ntheta-1) luminosity -= g*Rtht*dr*dphi;
                                    }
                                }
                            }
                        }
                    }
                } // izone
                
                if(mpi.getNumProcs() > 1) {
                    mpi.AllreduceSUM(shellAvgRho, shellAvgRho);
                    mpi.AllreduceSUM(shellAvgRho2, shellAvgRho2);
                    mpi.AllreduceSUM(shellAvgPress,shellAvgPress);
                    mpi.AllreduceSUM(shellAvgTgas, shellAvgTgas);
                    mpi.AllreduceSUM(shellAvgTgasZ0, shellAvgTgasZ0);
                    mpi.AllreduceSUM(shellAvgVphi, shellAvgVphi);
                    mpi.AllreduceSUM(shellAvgQ_MRI2, shellAvgQ_MRI2);
                    mpi.AllreduceSUM(shellAvgQ_MRI3, shellAvgQ_MRI3);
                    mpi.AllreduceSUM(shellAvgalpha1, shellAvgalpha1);
                    mpi.AllreduceSUM(shellMassFlux, shellMassFlux);
		            mpi.AllreduceSUM(shellAvgvzz,shellAvgvzz);
		            mpi.AllreduceSUM(shellAvgradV,shellAvgradV);

                    mpi.AllreduceSUM(shellAvgScaleHeight2, shellAvgScaleHeight2);
                    mpi.AllreduceSUM(shellAvgAlfvenSpeed, shellAvgAlfvenSpeed);
                    mpi.AllreduceSUM(shellAvgAlfvenSpeed2, shellAvgAlfvenSpeed2);
                    mpi.AllreduceSUM(shellAvgAlfvenSpeed3, shellAvgAlfvenSpeed3);
                    mpi.AllreduceSUM(shellAvgPtot0,shellAvgPtot0);
                    mpi.AllreduceSUM(shellAvgPmag0,shellAvgPmag0);
		            mpi.AllreduceSUM(shellAvgRho0,shellAvgRho0);
                    mpi.AllreduceSUM(faceArea4[irad],faceArea4[irad]);
                    for(int itheta = 0;itheta < ntheta; itheta++) {
                        mpi.AllreduceSUM(shellHeating[irad][itheta], shellHeating[irad][itheta]);
                        mpi.AllreduceSUM(shellStress[irad][itheta], shellStress[irad][itheta]);
                        mpi.AllreduceSUM(RadV[irad][itheta], RadV[irad][itheta]);
                        mpi.AllreduceSUM(Vz[irad][itheta], Vz[irad][itheta]);
                        mpi.AllreduceSUM(azAvgRho[itheta], azAvgRho[itheta]);
                        mpi.AllreduceSUM(azAvgv1[itheta], azAvgv1[itheta]);
                        mpi.AllreduceSUM(azAvgv3[itheta], azAvgv3[itheta]);
                        mpi.AllreduceSUM(azAvgb1[itheta], azAvgb1[itheta]);
                        mpi.AllreduceSUM(azAvgb3[itheta], azAvgb3[itheta]);
                    }
                    
                    if(isRadOn) {
                        mpi.AllreduceSUM(shellAvgPrad, shellAvgPrad );
                        mpi.AllreduceSUM(shellAvgPrad_Pgas, shellAvgPrad_Pgas);
                        mpi.AllreduceSUM(shellAvgTrad, shellAvgTrad);
                        mpi.AllreduceSUM(shellAvgTradZ0, shellAvgTradZ0);
                        mpi.AllreduceSUM(shellAvgErad,shellAvgErad);
                        for(int itheta = 0;itheta < ntheta; itheta++) {
                            mpi.AllreduceSUM(kappaTRho[irad][itheta], kappaTRho[irad][itheta]);
                            mpi.AllreduceSUM(kappaARho[irad][itheta], kappaARho[irad][itheta]);
                            mpi.AllreduceSUM(kappaScat[irad][itheta], kappaScat[irad][itheta]);
                            mpi.AllreduceSUM(shellFluxR[irad][itheta], shellFluxR[irad][itheta]);
                            mpi.AllreduceSUM(shellFluxZ[irad][itheta], shellFluxZ[irad][itheta]);
                            mpi.AllreduceSUM(azAvgErad[itheta], azAvgErad[itheta]);
                            mpi.AllreduceSUM(azAvgPratio[itheta], azAvgPratio[itheta]);
                            mpi.AllreduceSUM(azAvgvrad1[itheta], azAvgvrad1[itheta]);
                            mpi.AllreduceSUM(azAvgvrad3[itheta], azAvgvrad3[itheta]);
                        }
                    }
                }
                
                if(n>=nStartAverage) {
                    tsAvgQ_MRI2[irad]   += shellAvgQ_MRI2/shellAvgRho;
                    tsAvgQ_MRI3[irad]   += shellAvgQ_MRI3/shellAvgRho;
                    tsAvgalpha1[irad]   += shellAvgalpha1/shellAvgRho;
                    for(int itheta = 0;itheta < ntheta; itheta++) {
                        tsAvgRho[irad][itheta] += azAvgRho[itheta]/nphi;
                        tsAvgErad[irad][itheta] += azAvgErad[itheta]/nphi;
                        tsAvgPratio[irad][itheta] += azAvgPratio[itheta]/nphi;
                        tsAvgv1[irad][itheta] += azAvgv1[itheta]/nphi;
                        tsAvgv3[irad][itheta] += azAvgv3[itheta]/nphi;
                        tsAvgb1[irad][itheta] += azAvgb1[itheta]/nphi;
                        tsAvgb3[irad][itheta] += azAvgb3[itheta]/nphi;
                        tsAvgvrad1[irad][itheta] += azAvgvrad1[itheta]/nphi;
                        tsAvgvrad3[irad][itheta] += azAvgvrad3[itheta]/nphi;

                        tsAvgRadV[irad][itheta]  += RadV[irad][itheta]/nphi;
                        tsAvgVz[irad][itheta]    += Vz[irad][itheta]/nphi;
                    }
                    tAvgMassFlux[irad]  += shellMassFlux;
                }
                
                // Write the files for space-time diagrams
                if(rank == 0) {
                    char *append;
                    if( n == nStart && iRestart == 0 && irad == 0) {
                        append = "w";
                    } else {
                        append = "a";
                    }
                    double dr = (xmax-xmin)/nx;
                    if(ilogr > 0) dr *= rad[irad];
                    pp.writeDataFile("Sigma", append, 3, time,rad[irad], (munit/(lunit*lunit))*shellAvgRho/(boxz*rad[irad]*dr));
                    pp.writeDataFile("ScaleH", append, 3, time, rad[irad], sqrt(shellAvgScaleHeight2/shellAvgRho2));
                    pp.writeDataFile("Height", append, 3, time, rad[irad], rad[irad]*sqrt(shellAvgScaleHeight2/shellAvgRho2) );
                    pp.writeDataFile("Pressure",append, 3, time, rad[irad], shellAvgPress/shellAvgRho);
                    pp.writeDataFile("Tgas", append, 3, time, rad[irad], tempunit*shellAvgTgas/shellAvgRho2);
                    pp.writeDataFile("Vphi", append, 3, time, rad[irad], shellAvgVphi/shellAvgRho);
                    pp.writeDataFile("alpha", append, 3, time, rad[irad], shellAvgalpha1/shellAvgRho);
                    pp.writeDataFile("Ptot0", append, 3, time, rad[irad], shellAvgPtot0/nphi);
		            pp.writeDataFile("Rho0", append, 3, time, rad[irad], shellAvgRho0/nphi);
                    pp.writeDataFile("Mdot", append, 3, time, rad[irad], shellMassFlux);
		            pp.writeDataFile("vzz0", append, 3, time, rad[irad], shellAvgvzz);
		            pp.writeDataFile("radV0", append, 3, time, rad[irad], shellAvgradV);

		    
                    if(isMagFieldOn)
                        pp.writeDataFile("Pmag0", append, 3, time, rad[irad], shellAvgPmag0/nphi);
                    if(isRadOn) {
                        pp.writeDataFile("Prad", append, 3, time, rad[irad],shellAvgPrad/shellAvgRho);
                        pp.writeDataFile("Prad_Pgas", append, 3, time, rad[irad], shellAvgPrad_Pgas/shellAvgRho);
                        pp.writeDataFile("Trad", append, 3, time, rad[irad], tempunit*shellAvgTrad/shellAvgRho);
                        pp.writeDataFile("RadiationEnergy",append,3,time,rad[irad], shellAvgErad/shellAvgRho);
                    }
                }
                
                
            } //nrad
            mpi.AllreduceSUM(mass,mass);
            mpi.AllreduceSUM(luminosity,luminosity);
            
            for(int itheta = 0; itheta < ntheta; itheta++) {
                double zAvgRho        = 0.0;
                double zAvgBphi       = 0.0;
                double zAvgarad       = 0.0;
                double zAvgagas       = 0.0;
                double zAvgamag       = 0.0;
                double surfaceArea    = 0.0;
                
                for(int izone = 0; izone < rho.size(); izone++) {
                    if(fabs((posRTP[izone].getY() - theta[itheta])/theta[itheta]) < tol) {
                        if(posRTP[izone].getX() < radCut) {
                            double vol  = (*mesh).getZone(mzone[izone])->getVolume();
                            zAvgRho  += rho[izone]*sqrtg[izone]*vol;
                            surfaceArea += sqrtg[izone]*vol;
                        }
                        Vector3d dxCov = (*mesh).getZone(mzone[izone])->getVectorCovariantZoneLength();
                        if(fabs(posRTP[izone].getX() - radAccel) < 0.5*dxCov.getX()) {
                            //zAvgBphi += magfieldRPZ[izone].getY();
                            zAvgagas -= gradPgas[mzone[izone]].getY()/rho[izone];
                            zAvgamag -= gradPmag[mzone[izone]].getY()/rho[izone];
                            
                            double fourthirds = 4./3.;
                            double onethird   = 1./3.;
                            
                            Tensor4d metks,metinv;
                            metricks2.getMetric(posRTP[izone], metks);
                            metinv = metks.getInverse();
                            double g00 = metks.getTT();
                            double g01 = metks.getTX();
                            double g02 = metks.getTY();
                            double g03 = metks.getTZ();
                            double g11 = metks.getXX();
                            double g12 = metks.getXY();
                            double g13 = metks.getXZ();
                            double g22 = metks.getYY();
                            double g23 = metks.getYZ();
                            double g33 = metks.getZZ();
                            
                            double g00up = metinv.getTT();
                            double g01up = metinv.getTX();
                            double g02up = metinv.getTY();
                            double g03up = metinv.getTZ();
                            double g11up = metinv.getXX();
                            double g12up = metinv.getXY();
                            double g13up = metinv.getXZ();
                            double g22up = metinv.getYY();
                            double g23up = metinv.getYZ();
                            double g33up = metinv.getZZ();
                            
                            // radiation fields
                            double vrr   = radVelRTP[izone].getX();
                            double vrth  = radVelRTP[izone].getY();
                            double vrph  = radVelRTP[izone].getZ();
                            double ur0up = sqrt(-1./(metks(0,0)+2.*(metks(0,1)*vrr+metks(0,2)*vrth+metks(0,3)*vrph)+2.*(metks(1,2)*vrr*vrth+metks(2,3)*vrth*vrph+metks(1,3)*vrr*vrph)+metks(1,1)*vrr*vrr+metks(2,2)*vrth*vrth+metks(3,3)*vrph*vrph));
                            double ur1up = vrr*ur0up;
                            double ur2up = vrth*ur0up;
                            double ur3up = vrph*ur0up;
                            
                            double velr = velocityRTP[izone].getX();
                            double velt = velocityRTP[izone].getY();
                            double velph = velocityRTP[izone].getZ();
                            double u0up = boost[izone]/sqrtg[izone];
                            double u1up = velr*u0up;
                            double u2up = velt*u0up;
                            double u3up = velph*u0up;
                            double u0dn = g00*u0up + g01*u1up + g02*u2up + g03*u3up;
                            double u1dn = g01*u0up + g11*u1up + g12*u2up + g13*u3up;
                            double u2dn = g02*u0up + g12*u1up + g22*u2up + g23*u3up;
                            double u3dn = g03*u0up + g13*u1up + g23*u2up + g33*u3up;
                            
                            double Rmn00up = radiationEnergy[izone]*(fourthirds*ur0up*ur0up + onethird*g00up);
                            double Rmn01up = radiationEnergy[izone]*(fourthirds*ur0up*ur1up + onethird*g01up);
                            double Rmn02up = radiationEnergy[izone]*(fourthirds*ur0up*ur2up + onethird*g02up);
                            double Rmn03up = radiationEnergy[izone]*(fourthirds*ur0up*ur3up + onethird*g03up);
                            double Rmn11up = radiationEnergy[izone]*(fourthirds*ur1up*ur1up + onethird*g11up);
                            double Rmn12up = radiationEnergy[izone]*(fourthirds*ur1up*ur2up + onethird*g12up);
                            double Rmn13up = radiationEnergy[izone]*(fourthirds*ur1up*ur3up + onethird*g13up);
                            double Rmn22up = radiationEnergy[izone]*(fourthirds*ur2up*ur2up + onethird*g22up);
                            double Rmn23up = radiationEnergy[izone]*(fourthirds*ur2up*ur3up + onethird*g23up);
                            double Rmn33up = radiationEnergy[izone]*(fourthirds*ur3up*ur3up + onethird*g33up);
                            
                            double kapaP, kapaR, kaps, kap, arT4;
                            vector<double> opaca, opacs;
                            (*opacity).getOpacity("absorption", "planckMean", rho[izone], temperature[izone], opaca);
                            if(opaca.size() > 0) kapaP = opaca[0];
                            opaca.clear();
                            (*opacity).getOpacity("absorption", "rossMean", rho[izone], temperature[izone], opaca);
                            if(opaca.size() > 0) kapaR = opaca[0];
                            opaca.clear();
                            (*opacity).getOpacity("scattering", rho[izone], temperature[izone], opacs);
                            if(opacs.size() > 0) kaps  = opacs[0];
                            opacs.clear();
                            kap  = kapaR + kaps;
                            arT4 = aCon*pow(temperature[izone],4);
                            
                            double Rdoubledotu = (Rmn00up*u0dn + Rmn01up*u1dn + Rmn02up*u2dn + Rmn03up*u3dn)*u0dn + (Rmn01up*u0dn + Rmn11up*u1dn + Rmn12up*u2dn + Rmn13up*u3dn)*u1dn + (Rmn02up*u0dn + Rmn12up*u1dn + Rmn22up*u2dn + Rmn23up*u3dn)*u2dn + (Rmn03up*u0dn + Rmn13up*u1dn + Rmn23up*u2dn + Rmn33up*u3dn)*u3dn;
                            
                            double G0up = -rho[izone]*(kap*(Rmn00up*u0dn + Rmn01up*u1dn + Rmn02up*u2dn + Rmn03up*u3dn) + ((kap-kapaP)*Rdoubledotu + kapaP*arT4)*u0up);
                            double G1up = -rho[izone]*(kap*(Rmn01up*u0dn + Rmn11up*u1dn + Rmn12up*u2dn + Rmn13up*u3dn) + ((kap-kapaP)*Rdoubledotu + kapaP*arT4)*u1up);
                            double G2up = -rho[izone]*(kap*(Rmn02up*u0dn + Rmn12up*u1dn + Rmn22up*u2dn + Rmn23up*u3dn) + ((kap-kapaP)*Rdoubledotu + kapaP*arT4)*u2up);
                            double G3up = -rho[izone]*(kap*(Rmn03up*u0dn + Rmn13up*u1dn + Rmn23up*u2dn + Rmn33up*u3dn) + ((kap-kapaP)*Rdoubledotu + kapaP*arT4)*u3up);
                            if(isComptonOn) {
                                double trad = pow(radiationEnergy[izone]/aCon,0.25);
                                double coef = 4.*rho[izone]*kaps*kb*(temperature[izone]-trad)/(m_e*cee*cee)*Rdoubledotu;
                                G0up -= coef*u0up;
                                G1up -= coef*u1up;
                                G2up -= coef*u2up;
                                G3up -= coef*u3up;
                            }
                            double G3dn = g03*G0up + g13*G1up + g23*G2up + g33*G3up;
                            zAvgarad += G3dn/rho[izone];
                        }
                    }
                }
                if(mpi.getNumProcs() > 1) {
                    mpi.AllreduceSUM(zAvgRho, zAvgRho);
                    mpi.AllreduceSUM(zAvgBphi, zAvgBphi);
                    mpi.AllreduceSUM(zAvgarad, zAvgarad);
                    mpi.AllreduceSUM(zAvgagas, zAvgagas);
                    mpi.AllreduceSUM(zAvgamag, zAvgamag);
                    mpi.AllreduceSUM(surfaceArea, surfaceArea);
                }
                
                if(n >= nStartAverage) {
                    tzAvgRho[itheta] += zAvgRho/surfaceArea;
                    tzAvgarad[itheta] += zAvgarad/nphi;
                    tzAvgagas[itheta] += zAvgagas/nphi;
                    tzAvgamag[itheta] += zAvgamag/nphi;
                }
                
                if(rank == 0) {
                    char *append;
                    if( n == nStart && iRestart == 0 && itheta == 0) {
                        append = "w";
                    } else {
                        append = "a";
                    }
                    pp.writeDataFile("zAvgRho",append,2,time,zAvgRho/surfaceArea);
                    pp.writeDataFile("zAvgBphi",append,2,time,zAvgBphi/nphi);
                }
            } // nz
            
            if(rank == 0) {
                // Find photosphere and integrate quantities between top & bottom photospheres
                double rAvgPhotoZmax = 0., rAvgPhotoZmin = 0., nAvgRad = 0.;
                for(int irad =0;irad<nrad;irad++) {
                    char *append;
                    if( n == nStart && iRestart == 0 && irad == 0) {
                        append = "w";
                    } else {
                        append = "a";
                    }
                    double coolingRate = 0.0, heatingRate = 0.0, verticalStress = 0.0;
                    double tauT_b = 0.0, tauT_t = 0.0, tau_scat = 0.0, tauA_b = 0.0, tauA_t = 0.0, tau_Eff = 0.0;
                    double photoZmin, photoZmax,photoZmax_scat,photoZmax_abs,photoZmin_scat,photoZmin_abs, photoZmax_Eff;
                    for(int itheta = 0;itheta<ntheta/2;itheta++) {
                        tauT_t += kappaTRho[irad][itheta]/nphi;
                        
                        if(tauT_t < 1.){
                            tauT_t += kappaTRho[irad][itheta]/nphi;
                            photoZmax = rad[irad]*cos(theta[itheta]);
                        }
                        if(tau_scat < 1.) {
                            tau_scat      += kappaScat[irad][itheta]/nphi;
                            photoZmax_scat = rad[irad]*cos(theta[itheta]);
                        }
                        if(tauA_t < 1.){
                            tauA_t       += kappaARho[irad][itheta]/nphi;
                            photoZmax_abs = rad[irad]*cos(theta[itheta]);
                        }

                        if(tau_Eff < 1.){
                            tau_Eff      += sqrt(kappaARho[irad][itheta]*kappaScat[irad][itheta]/nphi/nphi);
                            photoZmax_Eff = rad[irad]*cos(theta[itheta]);
                        }
                        if((tauT_t > 1. || rad[irad]*cos(theta[itheta]) < photoZmax) && itheta > 0) {
                            verticalStress += shellStress[irad][itheta]/nphi;
                            heatingRate   += shellHeating[irad][itheta]/nphi;
                            coolingRate += (shellFluxZ[irad][itheta+1]-shellFluxZ[irad][itheta-1])/2./nphi;
                        }
                    }
                    if(n >= nStartAverage) {
                    tAvgPhotoZmax[irad] += photoZmax;
                    tAvgPhotoZmax_Scat[irad] += photoZmax_scat;
                    tAvgPhotoZmax_Abs[irad]  += photoZmax_abs;
                    tAvgPhotoZmax_Eff[irad]  += photoZmax_Eff;
                    }
                    for(int itheta=(ntheta-1);itheta>=ntheta/2;itheta--) {
                        cout << " " << rad[irad] << " " << rad[irad]*cos(theta[ntheta/2]) << endl ;

                        tauT_b += kappaTRho[irad][itheta]/nphi;
                        
                        if(tauT_b < 1.) {
                            tauT_t += kappaTRho[irad][itheta]/nphi;
                            photoZmin = rad[irad]*cos(theta[itheta]);
                        }
                        if(tau_scat < 1.){
                        tau_scat      += kappaScat[irad][itheta]/nphi;
                        photoZmin_scat = rad[irad]*cos(theta[itheta]);
                        }
                        if(tauA_b < 1.) {
                            tauA_b += kappaARho[irad][itheta]/nphi;
                            photoZmin_abs = rad[irad]*cos(theta[itheta]);
                        }
                        if((tauT_b > 1. || rad[irad]*cos(theta[itheta]) > photoZmin) && itheta < ntheta-1) {
                            verticalStress += shellStress[irad][itheta]/nphi;
                            heatingRate   += shellHeating[irad][itheta]/nphi;
                            coolingRate += (shellFluxZ[irad][itheta+1]-shellFluxZ[irad][itheta-1])/2./nphi;
                        }
                    }
                    if(rad[irad] < radCut) {
                        rAvgPhotoZmin += photoZmin;
                        rAvgPhotoZmax += photoZmax;
                        nAvgRad += 1.;
                    }
                    if(n >= nStartAverage) {
                        tAvgPhotoZmin[irad] += photoZmin;
                        tAvgPhotoZmin_Scat[irad] += photoZmin_scat;
                        tAvgPhotoZmin_Abs[irad] += photoZmin_abs;
                        tsAvgCoolingRate[irad] += coolingRate*(faceArea4[irad]/nphi);
                    }
                    pp.writeDataFile("CoolingRate",append,3,time,rad[irad],coolingRate);
                    pp.writeDataFile("HeatingRate",append,3,time,rad[irad],heatingRate);
                    pp.writeDataFile("VerticalStress",append,3,time,rad[irad],munit/(tunit*tunit)*verticalStress);
                } // irad
                
                char *append;
                if( n == nStart && iRestart == 0 ) {
                    append = "w";
                } else {
                    append = "a";
                }
                pp.writeDataFile("rAvgPhotoZmax",append,2,time,rAvgPhotoZmax/nAvgRad);
                pp.writeDataFile("rAvgPhotoZmin",append,2,time,rAvgPhotoZmin/nAvgRad);
                pp.writeDataFile("totalMass",append,2,time,mass);
                pp.writeDataFile("luminosity",append,2,time,luminosity);
            }
            
            timeold = time;
        } // end if
        
        
    } //time step
    
    /*write the time and shell averaged data*/
    if(rank == 0){
        for(int irad =0;irad<nrad;irad++) {
            char *append;
            if( irad == 0) {
                append = "w";
            } else {
                append = "a";
            }
            if(isRadOn) {
                pp.writeDataFile("tAvgPhotoZmax",append,2,rad[irad],tAvgPhotoZmax[irad]/nAverage);
                pp.writeDataFile("tAvgPhotoZmin",append,2,rad[irad],tAvgPhotoZmin[irad]/nAverage);
                pp.writeDataFile("tAvgPhotoZmax_Scat",append,2,rad[irad],tAvgPhotoZmax_Scat[irad]/nAverage);
                pp.writeDataFile("tAvgPhotoZmin_Scat",append,2,rad[irad],tAvgPhotoZmin_Scat[irad]/nAverage);
                pp.writeDataFile("tAvgPhotoZmax_Abs",append,2,rad[irad],tAvgPhotoZmax_Abs[irad]/nAverage);
                pp.writeDataFile("tAvgPhotoZmin_Abs",append,2,rad[irad],tAvgPhotoZmin_Abs[irad]/nAverage);
                pp.writeDataFile("tAvgPhotoZmax_Eff",append,2,rad[irad],tAvgPhotoZmax_Eff[irad]/nAverage);
                
            }
            if(isMagFieldOn)
            {
                pp.writeDataFile("MRIQ2",append,2,rad[irad],tsAvgQ_MRI2[irad]/nAverage);
                pp.writeDataFile("MRIQ3",append,2,rad[irad],tsAvgQ_MRI3[irad]/nAverage);
            }
            pp.writeDataFile("tsAvgalpha",append,2,rad[irad],tsAvgalpha1[irad]/nAverage);
            for(int itheta = 0;itheta < ntheta; itheta++) {
                if( irad == 0 && itheta == 0) {
                    append = "w";
                } else {
                    append = "a";
                }
                double rcostheta = rad[irad]*cos(theta[itheta]);
                double rsintheta = rad[irad]*sin(theta[itheta]);
                pp.writeDataFile("RadV_ting",append,3,rad[irad],rcostheta,tsAvgRadV[irad][itheta]);
                pp.writeDataFile("Vz_ting",append,3,rad[irad],rcostheta,tsAvgVz[irad][itheta]);

                pp.writeDataFile("tsAvgRho",append,3,rad[irad],rcostheta,tsAvgRho[irad][itheta]/nAverage);
                pp.writeDataFile("tsAvgErad",append,3,rad[irad],rcostheta,tsAvgErad[irad][itheta]/nAverage);
                pp.writeDataFile("tsAvgPratio",append,3,rad[irad],rcostheta,tsAvgPratio[irad][itheta]/nAverage);
                pp.writeDataFile("tsAvgStreamVel",append,4,rad[irad],rcostheta,tsAvgv1[irad][itheta]/nAverage,tsAvgv3[irad][itheta]/nAverage);
                if(isMagFieldOn)
                    pp.writeDataFile("tsAvgStreamMag",append,4,rad[irad],theta[itheta],tsAvgb1[irad][itheta]/nAverage,tsAvgb3[irad][itheta]/nAverage);
                pp.writeDataFile("tsAvgStreamRad",append,4,rad[irad],rcostheta,tsAvgvrad1[irad][itheta]/nAverage,tsAvgvrad3[irad][itheta]/nAverage);
            }
            
        }
        /*for(int itheta = 0;itheta <ntheta;itheta++) {
         char *append;
         if( itheta == 0) {
         append = "w";
         } else {
         append = "a";
         }
         pp.writeDataFile("tzAvgRho",append,2,theta[itheta],tzAvgRho[itheta]/nAverage);
         pp.writeDataFile("tzAvgarad",append,2,theta[itheta],tzAvgarad[itheta]/nAverage);
         pp.writeDataFile("tzAvgagas",append,2,theta[itheta],tzAvgagas[itheta]/nAverage);
         pp.writeDataFile("tzAvgamag",append,2,theta[itheta],tzAvgamag[itheta]/nAverage);
         }*/
    }
    mpi.finalize();
}

//Routine to calculate tensor components in the comoving frame.
void comovingTensor(Vector3d posCyl, double bhspin, Tensor4d gBL, Vector4d uCyl,
                    vector<vector<double> > TCyl, Tensor3d &Tout)
{
    vector<vector<double> > jacob(4), Tin(4), TinBL(4);
    Vector4d uKS = Vector4d(0.), uBL = Vector4d(0.), er = Vector4d(0.),
    etheta = Vector4d(0.), ephi = Vector4d(0.);
    for(int i = 0; i < 4; i++) {
        jacob[i].resize(4, 0.);
        Tin[i].resize(4, 0.);
        TinBL[i].resize(4, 0.);
    }
    
    double R = posCyl.getX();
    double z = posCyl.getZ();
    double rad = sqrt(R*R+z*z);
    
    //jacob[sph,cyl]
    jacob[0][0]=1.0;
    jacob[1][1]=R/rad;
    jacob[1][3]=z/rad;
    jacob[2][1]=z/rad/rad;
    jacob[2][3]=-R/rad/rad;
    jacob[3][2]=1.0;
    
    //loop to convert the Cyl 4-velocity to Sph 4-velocity
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            uKS[i] += jacob[i][j]*uCyl[j];
        }}
    
    //reset jacob variable
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            jacob[i][j] = 0.;
        }}
    
    //jacob[cyl,sph]
    jacob[0][0]=1.0;
    jacob[1][1]=R/rad;
    jacob[3][1]=z/rad;
    jacob[1][2]=z;
    jacob[3][2]=-R;
    jacob[2][3]=1.0;
    
    //loop to convert the Cyl Stress-Energy Tensor to Sph Stress Energy Tensor
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            for(int k = 0; k < 4; k++) {
                for(int l = 0; l < 4; l++) {
                    Tin[i][j] += jacob[k][i]*jacob[l][j]*TCyl[k][l];
                }}}}
    
    //reset jacob variable
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            jacob[i][j] = 0.;
        }}
    
    double Delta = rad*rad - 2.*rad + bhspin*bhspin;
    jacob[0][0]=1.0;
    jacob[0][1]=2.0*rad/Delta;
    jacob[1][1]=1.0;
    jacob[2][2]=1.0;
    jacob[3][1]=bhspin/Delta;
    jacob[3][3]=1.0;
    
    //loop to convert the KS Stress-Energy Tensor to BL Stress Energy Tensor
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            for(int k = 0; k < 4; k++) {
                for(int l = 0; l < 4; l++) {
                    TinBL[i][j] += jacob[k][i]*jacob[l][j]*Tin[k][l];
                }}}}
    
    //invert the jacobian matrix
    jacob[0][1] = -jacob[0][1];
    jacob[3][1] = -jacob[3][1];
    
    //loop to convert the KS 4-velocity to BL 4-velocity
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            uBL[i] += jacob[i][j]*uKS[j];
        }}
    
    //Tetrads (Beckwith, Hawley & Krolik 2008)
    double Kr=sqrt(abs(gBL(3,3)*uBL[3]*uBL[3] +
                       gBL(1,1)*uBL[1]*uBL[1] +
                       uBL[0]*(2.0*gBL(0,3)*uBL[3] + gBL(0,0)*uBL[0])));
    
    double Ktheta=sqrt(abs(gBL(3,3)*uBL[3]*uBL[3]+
                           gBL(1,1)*uBL[1]*uBL[1]+
                           gBL(2,2)*uBL[2]*uBL[2]+
                           uBL[0]*(2.0*gBL(0,3)*uBL[3]+gBL(0,0)*uBL[0])));
    
    double Kphi=sqrt(abs(gBL(3,3)*uBL[3]*uBL[3] +
                         uBL[0]*(2.0*gBL(0,3)*uBL[3] + gBL(0,0)*uBL[0])));
    
    er[0] = sqrt(gBL(1,1))*uBL[1]*uBL[0];
    er[1] = Kphi*Kphi/sqrt(gBL(1,1));
    er[3] = sqrt(gBL(1,1))*uBL[1]*uBL[3];
    er *= -1./(Kr*Kphi);
    
    etheta[0] = sqrt(gBL(2,2))*uBL[2]*uBL[0];
    etheta[1] = sqrt(gBL(2,2))*uBL[2]*uBL[1];
    etheta[2] = Kr*Kr/sqrt(gBL(2,2));
    etheta[3] = sqrt(gBL(2,2))*uBL[2]*uBL[3];
    etheta *= -1./(Kr*Ktheta);
    
    ephi[0] = gBL(3,3)*uBL[3] + gBL(0,3)*uBL[0];
    ephi[3] = -(gBL(0,3)*uBL[3] + gBL(0,0)*uBL[0]);
    ephi *= -1.0/(Kphi*sqrt(abs(-gBL(0,3)*gBL(0,3) + gBL(3,3)*gBL(0,0))));
    
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            Tout(0,2) += er[i]*ephi[j]*TinBL[i][j];
            Tout(0,1) += er[i]*etheta[j]*TinBL[i][j];
            Tout(1,2) += etheta[i]*ephi[j]*TinBL[i][j];
        }}
    
}


//End of main function
void setZoneConversion(vector<int> &mzone,
                       int  nx, int  ny, int  nz,
                       int  px, int  py, int  pz,
                       int opx, int opy, int opz)
{
    int lx = nx/px;
    int ly = ny/py;
    int lz = nz/pz;
    int pxr = opx/px;
    int pyr = opy/py;
    int pzr = opz/pz;
    int olx = lx/pxr;
    int oly = ly/pyr;
    int olz = lz/pzr;
    int ppZones = lx*ly*lz;
    int i = 0; int j = 0; int k = 0;
    for(int n = 0; n < ppZones; n++) {
        mzone[n] = k*lx*ly + j*lx + i;
        i+=1;
        if(i % olx == 0) {
            i-=olx;
            j+=1;
            if(j % oly == 0) {
                j-=oly;
                k+=1;
                if(k % olz == 0) {
                    k-=olz;
                    i+=olx;
                    if(i % lx == 0) {
                        i=0;
                        j+=oly;
                        if(j % ly == 0) {
                            j=0;
                            k+=olz;
                        }
                    }
                }
            }
        }
    }
}



