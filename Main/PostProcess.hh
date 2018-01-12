
#include <math.h>
#include <iostream>
#include <string.h>
#include <hdf5.h>
#include <vector>
#include <stdarg.h>
#include <cstdio>
#include <stdlib.h>     /* atoi */

using namespace std;
//========================================================================
// Cosmos++ post-processing class
//
// NOTE: geometry calculations (e.g., volume) assume Cartesian metric
//       and (x,y,z) ordered Cartesian Vectors (and nodepositions)
//       e.g., node 0 (x,y,z), node 1 (x,y,z), etc....
//========================================================================
class PPCLASS {
public:
    // Constructor/destructor
    PPCLASS(MPIClass mpi);
    ~PPCLASS() {};
    
    void readMasterFile(string rootDir, string fileStem = "Master");
    void readData(int dumpID);
    void readData(int dumpID, int px, int py, int pz,
                  int opx, int opy, int opz);
    void clearFields();
    
    void setZoneGeometry();
    void setBoxDimensions();
    void setVerbose(bool verbose);
    void setSphericalMetric();
    void setLogarithmicRadius(double r0, double eta0);
    void setNodePosition(vector<double> &nodepos);
    
    int  getScalarFieldIndex(string fieldName);
    int  getVectorFieldIndex(string fieldName);
    
    int  getParticleTypeIndex(string type);
    int  getParticleScalarFieldIndex(int n, string fieldName);
    int  getParticleVectorFieldIndex(int n, string fieldName);

    int  getDimension();
    int  getNumZones();
    int  getNumNodesInZone();
    int  getNumTimeSteps();
    
    vector<string> getScalarFieldNames();
    vector<string> getVectorFieldNames();
    
    vector<string> getParticleTypeNames();
    vector<string> getParticleScalarFieldNames(int n);
    vector<string> getParticleVectorFieldNames(int n);

    vector<double>  getGlobalBoxDimensions();
    vector<double>  getLocalBoxDimensions();
    
    vector<double>  getProblemDimensions();
    vector<double>  getZoneLengthPolar(int n);
    double dotProduct(vector<double>, vector<double>);
    
    double getZoneVolume(int n);
    void  getZoneVolume(vector<double> &vol);
    
    double getCurrentTime(int n);
    int   getCurrentCycle(int n);
    
    void  getDumpIDArray(vector<int>   &dmpid);
    void  getCycleArray(vector<int>    &cyc);
    void  getTimeArray(vector<double>   &tm);
    
    vector<double> getZoneLength(int n);
    void          getZoneLength(vector<double> &len);
    
    vector<double> getZoneCenter(int n);
    void          getZoneCenter(vector<double> &cen);
    
    int   getZoneID(vector<double> pos);
    int   getZoneID(vector<double> pos, int izone);
    
    void  getNeighborList(int izone, vector<int> &neighbors);
    
    void  getNodePosition(vector<double> &nodepos);
    void  getNodePosition(int izone, vector<double> &nodepos);
    
    void  getScalarField(string fieldname, vector<double> &field);
    void  getVectorField(string fieldname, vector<double> &field);
    
    void  getParticleScalarField(string type, string fieldname, vector<double> &field);
    void  getParticleVectorField(string type, string fieldname, vector<double> &field);

    void  writeDataFile(string fileName, char const *append, int n, ...);
    void  writeDataFileArray(string fileName, char const *append, int maxZone, double time, vector<double> data);
    
    double ran2(int idum);
    int   getRandomZoneID(int seed);
    
private:
    MPIClass mMPI;
    
    int mDimension;
    int mNumZones;
    int mNumNodesInZone;
    int mMetricType;
    int mNumDomains;
    int mNumScalarFields;
    int mNumVectorFields;
    int mNumParticleTypes;
    int mNumTimeSteps;
    
    string mRootDirName;
    
    bool mVerbose;
    bool mLogR;
    bool mGridDumpOn;
    
    double mR0;
    double mEta0;
    
    // Master file attributes
    // don't need to be resized before each read cycle
    
    vector<int>    mNumParticles;
    vector<int>    mNumParticleScalarFields;
    vector<int>    mNumParticleVectorFields;
    vector<int>    mDumpIDArray;
    vector<int>    mCycleArray;
    vector<double> mTimeArray;
    vector<string> mDataFileNameArray;
    vector<string> mGridFileNameArray;
    vector<string> mDirNameArray;
    
    vector<string> mScalarFieldNames;
    vector<string> mVectorFieldNames;
    
    vector<string> mParticleTypeNames;
    vector<vector<string> > mParticleScalarFieldNames;
    vector<vector<string> > mParticleVectorFieldNames;

    // data storage arrays
    // these must be resized/cleared before each read cycle
    
    vector<double>  mLocalBoxDimensions;
    vector<double>  mGlobalBoxDimensions;
    
    vector<double>  mZoneVolume;
    vector<double>  mZoneLength;
    vector<double>  mZoneCenter;
    vector<double>  mNodePosition;
    
    vector<vector<double> > mScalarFields;
    vector<vector<double> > mVectorFields;

    vector<vector<vector<double> > > mParticleScalarFields;
    vector<vector<vector<double> > > mParticleVectorFields;
};


//====================================================================
// CONSTRUCTOR
//====================================================================
PPCLASS::PPCLASS(MPIClass mpi) :
mMPI            (mpi),
mVerbose        (false),
mLogR           (false),
mGridDumpOn     (true),
mR0             (1.0),
mEta0           (1.0),
mNumTimeSteps   (0),
mNumZones       (0),
mNumDomains     (0),
mNumScalarFields(0),
mNumVectorFields(0),
mNumParticleTypes(0),
mNumParticles    (0),
mNumParticleScalarFields(0),
mNumParticleVectorFields(0),
mMetricType     (0)
{}

//====================================================================
// one-liner access functions
//====================================================================
int             PPCLASS::getDimension()           {return mDimension;}
int             PPCLASS::getNumZones()            {return mNumZones;}
int             PPCLASS::getNumTimeSteps()        {return mNumTimeSteps;}
int             PPCLASS::getNumNodesInZone()      {return mNumNodesInZone;}
int             PPCLASS::getCurrentCycle(int n)   {return mCycleArray[n];}

double           PPCLASS::getZoneVolume(int n)     {return mZoneVolume[n];}
double           PPCLASS::getCurrentTime(int n)    {return mTimeArray[n];}

void            PPCLASS::setVerbose(bool verbose) {mVerbose = verbose;}
void            PPCLASS::setSphericalMetric()     {mMetricType = 2;}

vector<string>  PPCLASS::getScalarFieldNames()    {return mScalarFieldNames;}
vector<string>  PPCLASS::getVectorFieldNames()    {return mVectorFieldNames;}

vector<string>  PPCLASS::getParticleTypeNames()           {return mParticleTypeNames;}
vector<string>  PPCLASS::getParticleScalarFieldNames(int n)    {return mParticleScalarFieldNames[n];}
vector<string>  PPCLASS::getParticleVectorFieldNames(int n)    {return mParticleVectorFieldNames[n];}

vector<double>   PPCLASS::getGlobalBoxDimensions() {return mGlobalBoxDimensions;}
vector<double>   PPCLASS::getLocalBoxDimensions()  {return mLocalBoxDimensions;}

void            PPCLASS::setLogarithmicRadius(double r0, double eta0) {
    mLogR = true;
    mR0   = r0;
    mEta0 = eta0;
}

//====================================================================
// access sized field functions
//====================================================================
void PPCLASS::getDumpIDArray(vector<int> &dmpid)
{
    dmpid.clear();
    dmpid.resize(mDumpIDArray.size());
    dmpid = mDumpIDArray;
}

void PPCLASS::getCycleArray(vector<int> &cyc)
{
    cyc.clear();
    cyc.resize(mCycleArray.size());
    cyc = mCycleArray;
}

void PPCLASS::getTimeArray(vector<double> &tm)
{
    tm.clear();
    tm.resize(mTimeArray.size());
    tm = mTimeArray;
}


void PPCLASS::getZoneLength(vector<double> &len)
{
    len.clear();
    len.resize(mDimension*mNumZones);
    len = mZoneLength;
}


vector<double> PPCLASS::getZoneLength(int izone)
{
    vector<double> len; //(mDimension);
    int offset = mDimension*izone;
    
    for(int i = 0; i < mDimension; i++) {
        len.push_back(mZoneLength[offset+i]);
    }
    return len;
}

vector<double>  PPCLASS::getProblemDimensions()
{
    double huge =  1.0e30;
    
    vector<double> dims;
    switch(mDimension)
    {
        case 2:
        {
            double min_r = huge, max_r = 0.0;
            double min_th = huge, max_th = 0.0;
            int min_r_index = 0, max_r_index =0, min_th_index = 0, max_th_index = 0;
            
            for(int n = 0; n < mNumZones; n++)
            {
                if(mZoneCenter[2*n] < min_r)
                {	min_r = mZoneCenter[2*n]; min_r_index = n; }
                if(mZoneCenter[2*n] > max_r)
                {	max_r = mZoneCenter[2*n]; max_r_index = n; }
                if(mZoneCenter[2*n+1] < min_th)
                {	min_th = mZoneCenter[2*n+1]; min_th_index = n; }
                if(mZoneCenter[2*n+1] > max_th)
                {	max_th = mZoneCenter[2*n+1]; max_th_index = n; }
            }
            
            vector<double> zone_len = this->getZoneLengthPolar(min_r_index);
            dims.push_back(min_r - zone_len[0]/2.0);
            
            zone_len = this->getZoneLengthPolar(max_r_index);
            dims.push_back(max_r + zone_len[0]/2.0);
            
            zone_len = this->getZoneLengthPolar(min_th_index);
            dims.push_back(min_th - zone_len[1]/2.0);
            
            zone_len = this->getZoneLengthPolar(max_th_index);
            dims.push_back(max_th + zone_len[1]/2.0);
        }
            
    }
    return dims;
}


vector<double> PPCLASS::getZoneLengthPolar(int n)
{
    vector<double> len;
    
    vector<double> x;
    vector<double> y;
    
    x.push_back(mNodePosition[8*n + 2]-mNodePosition[8*n + 0]);
    x.push_back(mNodePosition[8*n + 3]-mNodePosition[8*n + 1]);
    
    y.push_back(mNodePosition[8*n + 4]-mNodePosition[8*n + 6]);
    y.push_back(mNodePosition[8*n + 5]-mNodePosition[8*n + 7]);
    
    len.push_back(sqrt(dotProduct(x,x)));
    
    double cos_theta = dotProduct(x,y)/sqrt(dotProduct(x,x)*dotProduct(y,y));
    
    len.push_back(acos(cos_theta));
    
    return len;
}

double PPCLASS::dotProduct(vector<double> x, vector<double> y)
{
    int n = x.size();
    if(n != y.size())
        return 0;
    if(n < 1)
        return 0;
    
    double result = 0.0;
    for(int i = 0; i < n; i++)
        result += x[i]*y[i];
    return result;
}


void PPCLASS::getZoneCenter(vector<double> &cen)
{
    cen.clear();
    cen.resize(mDimension*mNumZones);
    cen = mZoneCenter;
}


vector<double> PPCLASS::getZoneCenter(int izone)
{
    vector<double> cen;
    int offset = mDimension*izone;
    
    for(int i = 0; i < mDimension; i++) {
        cen.push_back(mZoneCenter[offset+i]);
    }
    return cen;
}


void PPCLASS::getZoneVolume(vector<double> &vols)
{
    vols.clear();
    vols.resize(mNumZones);
    vols = mZoneVolume;
}


void PPCLASS::setNodePosition(vector<double> &nodepos)
{
    if(mDimension*mNumNodesInZone*mNumZones != nodepos.size()) {
        cout << "****ERROR: size of node list is incorrect" << endl;
        exit(0);
    }
    mNodePosition.clear();
    mNodePosition.resize(mDimension*mNumNodesInZone*mNumZones);
    mNodePosition = nodepos;
    
    // contruct zone geometry
    setZoneGeometry();
    setBoxDimensions();
}


void PPCLASS::getNodePosition(vector<double> &nodepos)
{
    if(mDimension*mNumNodesInZone*mNumZones != nodepos.size()) {
        cout << "****ERROR: size of node list is incorrect" << endl;
        exit(0);
    }
    nodepos.clear();
    nodepos.resize(mDimension*mNumNodesInZone*mNumZones);
    nodepos = mNodePosition;
}


void PPCLASS::getScalarField(string name, vector<double> &field)
{
    field.clear();
    field.resize(mNumZones);
    int index = getScalarFieldIndex(name);
    field     = mScalarFields[index];
}


void PPCLASS::getVectorField(string name, vector<double> &field)
{
    field.clear();
    field.resize(mDimension*mNumZones);
    int index = getVectorFieldIndex(name);
    field     = mVectorFields[index];
}

void PPCLASS::getParticleScalarField(string particleType, string fieldName, vector<double> &field)
{
    field.clear();
    int type = getParticleTypeIndex(particleType);
    field.resize(mNumParticles[type]);
    int index = getParticleScalarFieldIndex(type, fieldName);
    field     = mParticleScalarFields[type][index];
}


void PPCLASS::getParticleVectorField(string particleType, string fieldName, vector<double> &field)
{
    field.clear();
    int type = getParticleTypeIndex(particleType);
    field.resize(mDimension*mNumParticles[type]);
    int index = getParticleVectorFieldIndex(type, fieldName);
    field     = mParticleVectorFields[type][index];
}


//====================================================================
// clear all vector objects
//====================================================================
void PPCLASS::clearFields()
{
    mZoneLength.clear();
    mZoneCenter.clear();
    mZoneVolume.clear();
    mNodePosition.clear();
    
    mLocalBoxDimensions.clear();
    mGlobalBoxDimensions.clear();
    
    for(int n = 0; n < mScalarFields.size(); n++) mScalarFields[n].clear();
    for(int n = 0; n < mVectorFields.size(); n++) mVectorFields[n].clear();
    for(int m = 0; m < mParticleScalarFields.size(); m++) {
        for(int n = 0; n < mParticleScalarFields[m].size(); n++) mParticleScalarFields[m][n].clear();
        mParticleScalarFields[m].clear();
    }
    for(int m = 0; m < mParticleVectorFields.size(); m++) {
        for(int n = 0; n < mParticleVectorFields[m].size(); n++) mParticleVectorFields[m][n].clear();
        mParticleVectorFields[m].clear();
    }
    mScalarFields.clear();
    mVectorFields.clear();
    mParticleScalarFields.clear();
    mParticleVectorFields.clear();
}


//====================================================================
// retrieve index of field in vector lists
//
// key words follow the convention of Cosmos++, e.g.,
//     massDensity
//     energy
//     pressure
//     totalEnergy
//     boostFactor
//     velocity
//     momentum
//     magneticField
//====================================================================
int PPCLASS::getScalarFieldIndex(string fieldName)
{
    for(int n = 0; n < mScalarFieldNames.size(); n++) {
        if(fieldName == mScalarFieldNames[n]) return n;
    }
    
    if(mMPI.getProcID() == 0) {
        cout << endl;
        cout << "****WARNING: could not locate field index" << endl;
    }
    exit(0);
}

int PPCLASS::getVectorFieldIndex(string fieldName)
{
    for(int n = 0; n < mVectorFieldNames.size(); n++) {
        if(fieldName == mVectorFieldNames[n]) return n;
    }
    
    if(mMPI.getProcID() == 0) {
        cout << endl;
        cout << "****WARNING: could not locate field index" << endl;
    }
    exit(0);
}

int PPCLASS::getParticleTypeIndex(string typeName)
{
    for(int n = 0; n < mParticleTypeNames.size(); n++) {
        if(typeName == mParticleTypeNames[n]) return n;
    }
    
    if(mMPI.getProcID() == 0) {
        cout << endl;
        cout << "****WARNING: could not locate field index" << endl;
    }
    exit(0);
}

int PPCLASS::getParticleScalarFieldIndex(int type, string fieldName)
{
    for(int n = 0; n < mParticleScalarFieldNames[type].size(); n++) {
        if(fieldName == mParticleScalarFieldNames[type][n]) return n;
    }
    
    if(mMPI.getProcID() == 0) {
        cout << endl;
        cout << "****WARNING: could not locate field index" << endl;
    }
    exit(0);
}

int PPCLASS::getParticleVectorFieldIndex(int type, string fieldName)
{
    for(int n = 0; n < mParticleVectorFieldNames[type].size(); n++) {
        if(fieldName == mParticleVectorFieldNames[type][n]) return n;
    }
    
    if(mMPI.getProcID() == 0) {
        cout << endl;
        cout << "****WARNING: could not locate field index" << endl;
    }
    exit(0);
}


//====================================================================
// COMPUTE BOX DIMENSIONS
//====================================================================
void PPCLASS::setBoxDimensions()
{
    double huge =  1.0e30;
    
    double xmin  =  huge;
    double xmax  = -huge;
    double ymin  =  huge;
    double ymax  = -huge;
    double zmin  =  huge;
    double zmax  = -huge;
    double pi    = 2.0*asin(1.0);
    double delta = pi/16.0;
    
    int ncnt = 0;
    for(int izone = 0; izone < mNumZones; izone++) {
        
        if(mDimension == 3 && mMetricType == 2) { // 3D spherical only
            double r0,r1,h0,h1,p0,p1;
            double dr,dh,dp,ravg,havg,pavg;
            
            r0 = sqrt(mNodePosition[ncnt+ 0]*mNodePosition[ncnt+ 0] +
                      mNodePosition[ncnt+ 1]*mNodePosition[ncnt+ 1] +
                      mNodePosition[ncnt+ 2]*mNodePosition[ncnt+ 2]);
            r1 = sqrt(mNodePosition[ncnt+ 3]*mNodePosition[ncnt+ 3] +
                      mNodePosition[ncnt+ 4]*mNodePosition[ncnt+ 4] +
                      mNodePosition[ncnt+ 5]*mNodePosition[ncnt+ 5]);
            h0 = mNodePosition[ncnt+ 2]/r0;
            if(h0 >  1.0) h0 =  1.0;
            if(h0 < -1.0) h0 = -1.0;
            h0 = acos(h0);
            h1 = mNodePosition[ncnt+11]/r0;
            if(h1 >  1.0) h1 =  1.0;
            if(h1 < -1.0) h1 = -1.0;
            h1 = acos(h1);
            if(mNodePosition[ncnt+ 2] < 0.0) {
                p0 = (mNodePosition[ncnt+ 0]/
                      sqrt(mNodePosition[ncnt+ 0]*mNodePosition[ncnt+ 0] +
                           mNodePosition[ncnt+ 1]*mNodePosition[ncnt+ 1]));
                if(p0 >  1.0) p0 =  1.0;
                if(p0 < -1.0) p0 = -1.0;
                p0 = acos(p0);
                if(mNodePosition[ncnt+ 1] < 0.0) p0 = 2.0*pi - p0;
                p1 = (mNodePosition[ncnt+12]/
                      sqrt(mNodePosition[ncnt+12]*mNodePosition[ncnt+12] +
                           mNodePosition[ncnt+13]*mNodePosition[ncnt+13]));
                if(p1 >  1.0) p1 =  1.0;
                if(p1 < -1.0) p1 = -1.0;
                p1 = acos(p1);
                if(mNodePosition[ncnt+13] <= 0.0) p1 = 2.0*pi - p1;
            } else {
                p0 = (mNodePosition[ncnt+ 9]/
                      sqrt(mNodePosition[ncnt+ 9]*mNodePosition[ncnt+ 9] +
                           mNodePosition[ncnt+10]*mNodePosition[ncnt+10]));
                if(p0 >  1.0) p0 =  1.0;
                if(p0 < -1.0) p0 = -1.0;
                p0 = acos(p0);
                if(mNodePosition[ncnt+10] < 0.0) p0 = 2.0*pi - p0;
                p1 = (mNodePosition[ncnt+18]/
                      sqrt(mNodePosition[ncnt+18]*mNodePosition[ncnt+18] +
                           mNodePosition[ncnt+19]*mNodePosition[ncnt+19]));
                if(p1 >  1.0) p1 =  1.0;
                if(p1 < -1.0) p1 = -1.0;
                p1 = acos(p1);
                if(mNodePosition[ncnt+19] <= 0.0) p1 = 2.0*pi - p1;
            }
            
            if(r0 < xmin) xmin = r0;
            if(r1 > xmax) xmax = r1;
            if(h0 < ymin) ymin = h0;
            if(h1 > ymax) ymax = h1;
            if(p0 < zmin) zmin = p0;
            if(p1 > zmax) zmax = p1;
            
            ncnt += mDimension*mNumNodesInZone;
            
        } else {
            
            for(int n = 0; n < mNumNodesInZone; n++) {
                if(mNodePosition[ncnt] < xmin) xmin = mNodePosition[ncnt];
                if(mNodePosition[ncnt] > xmax) xmax = mNodePosition[ncnt];
                ncnt++;
                
                if(mDimension > 1) {
                    if(mNodePosition[ncnt] < ymin) ymin = mNodePosition[ncnt];
                    if(mNodePosition[ncnt] > ymax) ymax = mNodePosition[ncnt];
                    ncnt++;
                }
                
                if(mDimension > 2) {
                    if(mNodePosition[ncnt] > zmax) zmax = mNodePosition[ncnt];
                    if(mNodePosition[ncnt] < zmin) zmin = mNodePosition[ncnt];
                    ncnt++;
                }
            }
        }
    }
    
    mLocalBoxDimensions.push_back(xmin);
    mLocalBoxDimensions.push_back(xmax);
    
    if(mDimension > 1) {
        mLocalBoxDimensions.push_back(ymin);
        mLocalBoxDimensions.push_back(ymax);
    }
    
    if(mDimension > 2) {
        if(mMetricType == 2) {
            if((2.0*pi-zmax) < delta) zmax = 2.0*pi;
            if(zmin < delta) zmin = 0.0;
        }
        mLocalBoxDimensions.push_back(zmin);
        mLocalBoxDimensions.push_back(zmax);
    }
    
    if(mVerbose) {
        cout << "My ID: " << mMPI.getProcID() << " Local Box Dimensions (xmin, xmax, ymin, ymax, zmin, zmax)" << endl;
        cout << xmin << " " << xmax << endl;
        if(mDimension > 1) cout << ymin << " " << ymax << endl;
        if(mDimension > 2) cout << zmin << " " << zmax << endl;
    }
    
    if(mMPI.getNumProcs() > 1) {
        double inxmin = xmin;
        double inxmax = xmax;
        double inymin = ymin;
        double inymax = ymax;
        double inzmin = zmin;
        double inzmax = zmax;
        
        mMPI.AllreduceMIN(inxmin, xmin);
        mMPI.AllreduceMAX(inxmax, xmax);
        
        if(mDimension > 1) {
            mMPI.AllreduceMIN(inymin, ymin);
            mMPI.AllreduceMAX(inymax, ymax);
        }
        
        if(mDimension > 2) {
            mMPI.AllreduceMIN(inzmin, zmin);
            mMPI.AllreduceMAX(inzmax, zmax);
        }
    }
    
    mGlobalBoxDimensions.push_back(xmin);
    mGlobalBoxDimensions.push_back(xmax);
    
    if(mDimension > 1) {
        mGlobalBoxDimensions.push_back(ymin);
        mGlobalBoxDimensions.push_back(ymax);
    }
    
    if(mDimension > 2) {
        mGlobalBoxDimensions.push_back(zmin);
        mGlobalBoxDimensions.push_back(zmax);
    }
    
}


//====================================================================
// RETURN ZONE ID CONTAINING PARTICLE POSITION
//====================================================================
int PPCLASS::getZoneID(vector<double> pos)
{
    double huge     = 1.0e30;
    double pi       = 2.0*asin(1.0);
    int   zoneID   = -1;
    bool  isInside = true;
    
    // first check if position is within local box bounds
    int n = 0;
    while( n < mLocalBoxDimensions.size() ) {
        if(pos[n+0] < mLocalBoxDimensions[n+0]) isInside = false;
        if(pos[n+1] > mLocalBoxDimensions[n+1]) isInside = false;
        n += 2;
    }
    
    if(!isInside) {
        return -1;
    } else {
        isInside = false;
    }
    
    int  ncnt  = 0;
    int  izone = 0;
    while(izone < mNumZones && !isInside) {
        double xmin =  huge;
        double xmax = -huge;
        double ymin =  huge;
        double ymax = -huge;
        double zmin =  huge;
        double zmax = -huge;
        
        for(int n = 0; n < mNumNodesInZone; n++) {
            if(mNodePosition[ncnt] < xmin) xmin = mNodePosition[ncnt];
            if(mNodePosition[ncnt] > xmax) xmax = mNodePosition[ncnt];
            ncnt++;
            
            if(mDimension > 1) {
                if(mNodePosition[ncnt] < ymin) ymin = mNodePosition[ncnt];
                if(mNodePosition[ncnt] > ymax) ymax = mNodePosition[ncnt];
                ncnt++;
            }
            
            if(mDimension > 2) {
                if(mNodePosition[ncnt] < zmin) zmin = mNodePosition[ncnt];
                if(mNodePosition[ncnt] > zmax) zmax = mNodePosition[ncnt];
                ncnt++;
            }
        }
        
        isInside = true;
        if(pos[0] < xmin || pos[1] > xmax) isInside = false;
        
        if(mDimension > 1) {
            if(pos[2] < ymin || pos[3] > ymax) isInside = false;
        }
        
        if(mDimension > 2) {
            if(pos[4] < zmin || pos[5] > zmax) isInside = false;
        }
        
        if(isInside) zoneID = izone;
        izone++;
    }
    if(!isInside) {
        cout << pos[0] << " " << pos[2] << " " << pos[4] << endl;
        cout << mLocalBoxDimensions[0] << " " << mLocalBoxDimensions[2] << " " << mLocalBoxDimensions[4] << endl;
    }
    
    return zoneID;
}



//====================================================================
// RETURN ZONE ID CONTAINING PARTICLE POSITION in NEIGHBORHOOD OF ZONE
//====================================================================
int PPCLASS::getZoneID(vector<double> pos, int myZoneID)
{
    double huge     = 1.0e30;
    double pi       = 2.0*asin(1.0);
    int   zoneID   = -1;
    bool  isInside = true;
    
    // first check if position is within local box bounds
    int n = 0;
    while( n < mLocalBoxDimensions.size() ) {
        if(pos[n+0] < mLocalBoxDimensions[n+0]) isInside = false;
        if(pos[n+1] > mLocalBoxDimensions[n+1]) isInside = false;
        n += 2;
    }
    
    if(!isInside) {
        return -1;
    } else {
        isInside = false;
    }
    
    
    double xmin =  huge;
    double xmax = -huge;
    double ymin =  huge;
    double ymax = -huge;
    double zmin =  huge;
    double zmax = -huge;
    
    // next check current zone
    int izone  = myZoneID;
    int offset = mNumNodesInZone*mDimension*izone;
    
    int ncnt = offset;
    for(int n = 0; n < mNumNodesInZone; n++) {
        if(mNodePosition[ncnt] < xmin) xmin = mNodePosition[ncnt];
        if(mNodePosition[ncnt] > xmax) xmax = mNodePosition[ncnt];
        ncnt++;
        
        if(mDimension > 1) {
            if(mNodePosition[ncnt] < ymin) ymin = mNodePosition[ncnt];
            if(mNodePosition[ncnt] > ymax) ymax = mNodePosition[ncnt];
            ncnt++;
        }
        
        if(mDimension > 2) {
            if(mNodePosition[ncnt] < zmin) zmin = mNodePosition[ncnt];
            if(mNodePosition[ncnt] > zmax) zmax = mNodePosition[ncnt];
            ncnt++;
            if(mMetricType == 2) {
                if(zmax-zmin > pi) {
                    if(mNodePosition[ncnt] > pi) {
                        double zmaxt = zmax;
                        zmax = max(zmaxt, mNodePosition[ncnt]);
                        zmin = min(zmaxt, mNodePosition[ncnt]);
                    }
                    if(mNodePosition[ncnt] > pi) {
                        double zmint = zmin;
                        zmax = max(zmint, mNodePosition[ncnt]);
                        zmin = min(zmint, mNodePosition[ncnt]);
                    }
                }
            }
            
        }
    }
    
    isInside = true;
    if(pos[0] < xmin || pos[1] > xmax) isInside = false;
    
    if(mDimension > 1) {
        if(pos[2] < ymin || pos[3] > ymax) isInside = false;
    }
    
    if(mDimension > 2) {
        if(pos[4] < zmin || pos[5] > zmax) isInside = false;
    }
    
    if(isInside) return izone;
    
    
    // finally check neighbor zones
    vector<int> neighbors;
    getNeighborList(myZoneID, neighbors);
    
    int ii = 0;
    while(ii < neighbors.size() && !isInside) {
        int izone  = neighbors[ii];
        int offset = mNumNodesInZone*mDimension*izone;
        
        int ncnt = offset;
        for(int n = 0; n < mNumNodesInZone; n++) {
            if(mNodePosition[ncnt] < xmin) xmin = mNodePosition[ncnt];
            if(mNodePosition[ncnt] > xmax) xmax = mNodePosition[ncnt];
            ncnt++;
            
            if(mDimension > 1) {
                if(mNodePosition[ncnt] < ymin) ymin = mNodePosition[ncnt];
                if(mNodePosition[ncnt] > ymax) ymax = mNodePosition[ncnt];
                ncnt++;
            }
            
            if(mDimension > 2) {
                if(mNodePosition[ncnt] < zmin) zmin = mNodePosition[ncnt];
                if(mNodePosition[ncnt] > zmax) zmax = mNodePosition[ncnt];
                ncnt++;
                if(mMetricType == 2) {
                    if(zmax-zmin > pi) {
                        if(mNodePosition[ncnt] > pi) {
                            double zmaxt = zmax;
                            zmax = max(zmaxt, mNodePosition[ncnt]);
                            zmin = min(zmaxt, mNodePosition[ncnt]);
                        }
                        if(mNodePosition[ncnt] > pi) {
                            double zmint = zmin;
                            zmax = max(zmint, mNodePosition[ncnt]);
                            zmin = min(zmint, mNodePosition[ncnt]);
                        }
                    }
                }
                
            }
        }
        
        isInside = true;
        if(pos[0] < xmin || pos[1] > xmax) isInside = false;
        
        if(mDimension > 1) {
            if(pos[2] < ymin || pos[3] > ymax) isInside = false;
        }
        
        if(mDimension > 2) {
            if(pos[4] < zmin || pos[5] > zmax) isInside = false;
        }
        
        if(isInside) zoneID = izone;
        ii++;
    }
    
    return zoneID;
}




//====================================================================
// COMPUTE CELL VOLUMES, LENGTHS
//====================================================================
void PPCLASS::setZoneGeometry()
{
    double dx,dy,dz,xavg,yavg,zavg;
    double pi = 2.0*asin(1.0);
    mZoneLength.resize( mScalarFields[0].size()*mDimension );
    mZoneCenter.resize( mScalarFields[0].size()*mDimension );
    mZoneVolume.resize( mScalarFields[0].size() );
    
    if(mDimension == 3)
    {
        int ncnt = 0;
        for(int n = 0; n < mNumZones; n++) {
            if(mMetricType == 2) {
                double r0,r1,h0,h1,p0,p1;
                double dr,dh,dp,ravg,havg,pavg;
                
                r0 = sqrt(mNodePosition[ncnt+ 0]*mNodePosition[ncnt+ 0] +
                          mNodePosition[ncnt+ 1]*mNodePosition[ncnt+ 1] +
                          mNodePosition[ncnt+ 2]*mNodePosition[ncnt+ 2]);
                r1 = sqrt(mNodePosition[ncnt+ 3]*mNodePosition[ncnt+ 3] +
                          mNodePosition[ncnt+ 4]*mNodePosition[ncnt+ 4] +
                          mNodePosition[ncnt+ 5]*mNodePosition[ncnt+ 5]);
                h0 = mNodePosition[ncnt+ 2]/r0;
                if(h0 >  1.0) h0 =  1.0;
                if(h0 < -1.0) h0 = -1.0;
                h0 = acos(h0);
                h1 = mNodePosition[ncnt+11]/r0;
                if(h1 >  1.0) h1 =  1.0;
                if(h1 < -1.0) h1 = -1.0;
                h1 = acos(h1);
                if(mNodePosition[ncnt+ 2] < 0.0) {
                    p0 = (mNodePosition[ncnt+ 0]/
                          sqrt(mNodePosition[ncnt+ 0]*mNodePosition[ncnt+ 0] +
                               mNodePosition[ncnt+ 1]*mNodePosition[ncnt+ 1]));
                    if(p0 >  1.0) p0 =  1.0;
                    if(p0 < -1.0) p0 = -1.0;
                    p0 = acos(p0);
                    if(mNodePosition[ncnt+ 1] < 0.0) p0 = 2.0*pi - p0;
                    p1 = (mNodePosition[ncnt+12]/
                          sqrt(mNodePosition[ncnt+12]*mNodePosition[ncnt+12] +
                               mNodePosition[ncnt+13]*mNodePosition[ncnt+13]));
                    if(p1 >  1.0) p1 =  1.0;
                    if(p1 < -1.0) p1 = -1.0;
                    p1 = acos(p1);
                    if(mNodePosition[ncnt+13] <= 0.0) p1 = 2.0*pi - p1;
                } else {
                    p0 = (mNodePosition[ncnt+ 9]/
                          sqrt(mNodePosition[ncnt+ 9]*mNodePosition[ncnt+ 9] +
                               mNodePosition[ncnt+10]*mNodePosition[ncnt+10]));
                    if(p0 >  1.0) p0 =  1.0;
                    if(p0 < -1.0) p0 = -1.0;
                    p0 = acos(p0);
                    if(mNodePosition[ncnt+10] < 0.0) p0 = 2.0*pi - p0;
                    p1 = (mNodePosition[ncnt+18]/
                          sqrt(mNodePosition[ncnt+18]*mNodePosition[ncnt+18] +
                               mNodePosition[ncnt+19]*mNodePosition[ncnt+19]));
                    if(p1 >  1.0) p1 =  1.0;
                    if(p1 < -1.0) p1 = -1.0;
                    p1 = acos(p1);
                    if(mNodePosition[ncnt+19] <= 0.0) p1 = 2.0*pi - p1;
                }
                
                ravg = 0.5*(r0+r1);
                havg = 0.5*(h0+h1);
                pavg = 0.5*(p0+p1);
                dr = r1-r0;
                dh = (h1-h0);
                dp = (p1-p0);
                //cout << ravg << " " << havg << " " << pavg << endl;
                
                mZoneLength[n*mDimension+0] = dr;
                mZoneLength[n*mDimension+1] = dh;
                mZoneLength[n*mDimension+2] = dp;
                
                mZoneVolume[n] = dr*dh*dp;
                
                mZoneCenter[n*mDimension+0] = ravg;
                mZoneCenter[n*mDimension+1] = havg;
                mZoneCenter[n*mDimension+2] = pavg;
                if(mLogR) {
                    mZoneVolume[n] /= ravg;
                    mZoneLength[n*mDimension+0] /= ravg;
                }
                
            } else {
                xavg = 0.125*(  mNodePosition[ncnt+ 0] + mNodePosition[ncnt+ 3]
                              + mNodePosition[ncnt+ 6] + mNodePosition[ncnt+ 9]
                              + mNodePosition[ncnt+12] + mNodePosition[ncnt+15]
                              + mNodePosition[ncnt+18] + mNodePosition[ncnt+21] );
                
                yavg = 0.125*(  mNodePosition[ncnt+ 1] + mNodePosition[ncnt+ 4]
                              + mNodePosition[ncnt+ 7] + mNodePosition[ncnt+10]
                              + mNodePosition[ncnt+13] + mNodePosition[ncnt+16]
                              + mNodePosition[ncnt+19] + mNodePosition[ncnt+22] );
                
                zavg = 0.125*(  mNodePosition[ncnt+ 2] + mNodePosition[ncnt+ 5]
                              + mNodePosition[ncnt+ 8] + mNodePosition[ncnt+11]
                              + mNodePosition[ncnt+14] + mNodePosition[ncnt+17]
                              + mNodePosition[ncnt+20] + mNodePosition[ncnt+23] );
                
                dx = 0.25*(  fabs(mNodePosition[ncnt+ 0] - xavg)
                           + fabs(mNodePosition[ncnt+ 3] - xavg)
                           + fabs(mNodePosition[ncnt+ 6] - xavg)
                           + fabs(mNodePosition[ncnt+ 9] - xavg)
                           + fabs(mNodePosition[ncnt+12] - xavg)
                           + fabs(mNodePosition[ncnt+15] - xavg)
                           + fabs(mNodePosition[ncnt+18] - xavg)
                           + fabs(mNodePosition[ncnt+21] - xavg) );
                
                dy = 0.25*(  fabs(mNodePosition[ncnt+ 1] - yavg)
                           + fabs(mNodePosition[ncnt+ 4] - yavg)
                           + fabs(mNodePosition[ncnt+ 7] - yavg)
                           + fabs(mNodePosition[ncnt+10] - yavg)
                           + fabs(mNodePosition[ncnt+13] - yavg)
                           + fabs(mNodePosition[ncnt+16] - yavg)
                           + fabs(mNodePosition[ncnt+19] - yavg)
                           + fabs(mNodePosition[ncnt+22] - yavg) );
                
                dz = 0.25*(  fabs(mNodePosition[ncnt+ 2] - zavg)
                           + fabs(mNodePosition[ncnt+ 5] - zavg)
                           + fabs(mNodePosition[ncnt+ 8] - zavg)
                           + fabs(mNodePosition[ncnt+11] - zavg)
                           + fabs(mNodePosition[ncnt+14] - zavg)
                           + fabs(mNodePosition[ncnt+17] - zavg)
                           + fabs(mNodePosition[ncnt+20] - zavg)
                           + fabs(mNodePosition[ncnt+23] - zavg) );
                
                mZoneVolume[n] = dx*dy*dz;
                
                mZoneLength[n*mDimension+0] = dx;
                mZoneLength[n*mDimension+1] = dy;
                mZoneLength[n*mDimension+2] = dz;
                
                mZoneCenter[n*mDimension+0] = xavg;
                mZoneCenter[n*mDimension+1] = yavg;
                mZoneCenter[n*mDimension+2] = zavg;
            }
            
            ncnt += mDimension*mNumNodesInZone;
        }
        
    } else if(mDimension == 2) {
        
        int ncnt = 0;
        for(int n = 0; n < mNumZones; n++) {
            xavg = 0.25*(  mNodePosition[ncnt+0] + mNodePosition[ncnt+2]
                         + mNodePosition[ncnt+4] + mNodePosition[ncnt+6] );
            
            yavg = 0.25*(  mNodePosition[ncnt+1] + mNodePosition[ncnt+3]
                         + mNodePosition[ncnt+5] + mNodePosition[ncnt+7] );
            
            dx = 0.5*(  fabs(mNodePosition[ncnt+0] - xavg)
                      + fabs(mNodePosition[ncnt+2] - xavg)
                      + fabs(mNodePosition[ncnt+4] - xavg)
                      + fabs(mNodePosition[ncnt+6] - xavg) );
            
            dy = 0.5*(  fabs(mNodePosition[ncnt+1] - yavg)
                      + fabs(mNodePosition[ncnt+3] - yavg)
                      + fabs(mNodePosition[ncnt+5] - yavg)
                      + fabs(mNodePosition[ncnt+7] - yavg) );
            
            mZoneVolume[n] = dx*dy;
            
            mZoneLength[n*mDimension+0] = dx;
            mZoneLength[n*mDimension+1] = dy;
            
            
            if(mMetricType == 2) {				// polar coordinates
                double r = sqrt(xavg*xavg + yavg*yavg);
                double theta = pi/2 - atan2(yavg, xavg);
                // the reason we take pi/2 less than this is that we want theta measured over from y-axis, since the jet will be
                // symmetric about the y-axis
                
                mZoneCenter[n*mDimension+0] = r;
                mZoneCenter[n*mDimension+1] = theta;
                
            }
            else {
                mZoneCenter[n*mDimension+0] = xavg;
                mZoneCenter[n*mDimension+1] = yavg;
            }
            
            ncnt += mDimension*mNumNodesInZone;
        }
        
    } else if(mDimension == 1) {
        
        int ncnt = 0;
        for(int n = 0; n < mNumZones; n++) {
            xavg = 0.5*(mNodePosition[ncnt+0] + mNodePosition[ncnt+1]);
            
            dx = fabs(mNodePosition[ncnt+0] - xavg)
            + fabs(mNodePosition[ncnt+1] - xavg);
            
            mZoneVolume[n] = dx;
            mZoneLength[n] = dx;
            mZoneCenter[n] = xavg;
            
            ncnt += mDimension*mNumNodesInZone;
        }
        
    } else {
        if(mMPI.getProcID() == 0) {
            cout << "****ERROR-volume: mDimension is incorrect" << endl;
        }
        exit(0);
    }
}

//====================================================================
// Generate a psuedo-random number with a period of 10^18 using the
// ran2 algorithm from Numerical Methods
// RNMX should approximate largest doubleing value less than 1
//====================================================================

double PPCLASS::ran2(int idum)
{
    const int IM1=2147483563, IM2=2147483399;
    const int IA1=40014, IA2=40692, IQ1=53668, IQ2=52774;
    const int IR1=12211, IR2=3791, NTAB=32, IMM1=IM1 - 1;
    const int NDIV = 1 + IMM1/NTAB;
    const double EPS=3.0e-16, RNMX=1.0 - EPS, AM=1.0/double(IM1);
    static int idum2=123456789, iy=0;
    static vector<int> iv(NTAB);
    int j,k;
    double temp;
    
    if(idum <= 0) {                    //Initialize
        idum=(idum==0 ? 1 : -idum);      //Be sure to prevent idum=0
        idum2=idum;
        for (j=NTAB+7;j>=0;j--) {        //Load the shuffle table (after 8 warmups)
            k=idum/IQ1;
            idum=IA1*(idum-k*IQ1)-k*IR1;
            if (idum < 0) idum += IM1;
            if (j < NTAB) iv[j]=idum;
        }
        iy=iv[0];
    }
    k=idum/IQ1;                        //Start here when not initializing.
    idum=IA1*(idum-k*IQ1)-k*IR1;       //Compute idum=(IA1*idum)%IM1 w/o overflow
    if (idum < 0) idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;     //Compute idum2=(IA2*idum)%IM2 likewise
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;                         //Will be in range 0..NTAB-1
    iy=iv[j] - idum2;                  //Here idum is shuffled, idum and idum2
    iv[j] = idum;                      //are combined to generate output
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;  //because users don't expect endpoints
    else return temp;
}


//====================================================================
// READ MASTER FILE
//====================================================================
void PPCLASS::readMasterFile(string rootDir, string fileStem)
{
    FILE *fp;
    int ndims,ndomains,nfields,dummy,timesteps;
    char descriptor[100],fieldname[30],dirname[30];
    
    mRootDirName    = rootDir;
    //string fileName = rootDir + string("out-Master.cosmos++");
    string fileName = rootDir + string("out-") + fileStem + string(".cosmos++");
    
    if((fp = fopen(fileName.data(), "r")) == NULL) {
        cout<<"Error: Could not open master file "<<fileName<<endl;
        exit(1);
    }
    
    if(mMPI.getProcID() == 0 && mVerbose) {
        cout << "\nReading Master file " << fileName.data() << endl;
    }
    
    // number of domains/directories
    fscanf(fp, "%s %i", &descriptor, &ndomains);
    mNumDomains = ndomains;
    if(mMPI.getProcID() == 0 && mVerbose) printf("%s %i\n", &descriptor, ndomains);
    
    // directory names
    fscanf(fp, "%s", &descriptor);
    if(mMPI.getProcID() == 0 && mVerbose) printf("%s\n", &descriptor);
    for(int i = 0; i < ndomains; i++) {
        fscanf(fp, "%s", &dirname);
        mDirNameArray.push_back(dirname);
        if(mMPI.getProcID() == 0 && mVerbose) printf("%s\n", &dirname);
    }
    
    // number of dimensions
    fscanf(fp, "%s %i", &descriptor, &ndims);
    mDimension = ndims;
    if(mMPI.getProcID() == 0 && mVerbose) printf("%s %i\n", &descriptor, ndims);
    
    mNumNodesInZone = 8;
    if(ndims == 1) mNumNodesInZone = 2;
    if(ndims == 2) mNumNodesInZone = 4;
    
    // number of output fields
    fscanf(fp, "%s %i", &descriptor, &nfields);
    if(mMPI.getProcID() == 0 && mVerbose) printf("%s %i\n", &descriptor, nfields);
    
    // field names
    mNumScalarFields = 0;
    mNumVectorFields = 0;
    for(int i = 0; i < nfields; i++) {
        fscanf(fp, "%s %s", &descriptor, &fieldname);
        
        if(string(descriptor) == string("#SCALARFIELD")) {
            mNumScalarFields++;
            mScalarFieldNames.push_back( string(fieldname) );
            
        } else if(string(descriptor) == string("#VECTORFIELD")) {
            if(string(fieldname) != string("nodePosition")) {
                mNumVectorFields++;
                mVectorFieldNames.push_back( string(fieldname) );
            }
            
        } else {
            if(mMPI.getProcID() == 0) {
                cout << "****WARNING: descriptor "
                << string(descriptor) << " is not defined" << endl;
            }
            exit(0);
        }
        
        if(mMPI.getProcID() == 0 && mVerbose) printf("%s %s\n", &descriptor, &fieldname);
    }

    fscanf(fp, "%s %i", &descriptor, &dummy);
    if(string(descriptor) == string("#NUMBEROFPARTICLETYPES")) {
        // number of particles
        mNumParticleTypes = dummy;
        if(mMPI.getProcID() == 0 && mVerbose) printf("%s %i\n", &descriptor, mNumParticleTypes);
        mNumParticles.resize(mNumParticleTypes);
        mNumParticleScalarFields.resize(mNumParticleTypes);
        mNumParticleVectorFields.resize(mNumParticleTypes);
        mParticleScalarFields.resize(mNumParticleTypes);
        mParticleVectorFields.resize(mNumParticleTypes);
        mParticleScalarFieldNames.resize(mNumParticleTypes);
        mParticleVectorFieldNames.resize(mNumParticleTypes);
        for(int type = 0; type < mNumParticleTypes; type++) {
            fscanf(fp, "%s %s", &descriptor, &fieldname);
            mParticleTypeNames.push_back( string(fieldname) );
                
            for(int j = 0; j < 50; j++) {
                fscanf(fp, "%s %s", &descriptor, &fieldname);

                if(string(descriptor) == string("#PARTICLESCALARFIELD")) {
                    mNumParticleScalarFields[type]++;
                    mParticleScalarFieldNames[type].push_back( string(fieldname) );
                    
                } else if(string(descriptor) == string("#PARTICLEVECTORFIELD")) {
                    mNumParticleVectorFields[type]++;
                    mParticleVectorFieldNames[type].push_back( string(fieldname) );
                    
                } else {
                    timesteps = atoi( fieldname );
                    break;
                }
                
                if(mMPI.getProcID() == 0 && mVerbose) printf("%s %s\n", &descriptor, &fieldname);
            }
            mParticleScalarFields[type].resize(mNumParticleScalarFields[type]);
            mParticleVectorFields[type].resize(mNumParticleVectorFields[type]);
        }
    } else {
        timesteps = dummy;
    }
    
    
    // descriptor for timestep count
    if(mMPI.getProcID() == 0 && mVerbose) printf("%s %i\n", &descriptor, timesteps);
    mNumTimeSteps = timesteps;
    
    // descriptor for sequential output
    fscanf(fp, "%s\n", &descriptor); // the \n is important here because you must complete this line before testing the length of the next line (below)
    if(mMPI.getProcID() == 0 && mVerbose) printf("%s\n", &descriptor);
    
    // read in dump number, cycle, time, number of total zones,
    // number of internal zones, and dump file names, in that order
    char  fname[30], gname[30];
    float time;
    int   dump, cycle;
    int   scheck = 0;
    
    // Finding out if grid data is dumped separately from other data by
    // looking at how many objects are read by fscanf.  Need to reset the
    // file stream after this to properly read information.
    fpos_t position_before;
    fgetpos(fp, &position_before);
    char linestring[100];
    fgets(linestring, 100, fp);
    if(sscanf(linestring, "%i %i %e %s %s", &dump, &cycle, &time, &fname, &gname) != 5) mGridDumpOn = false;
    fsetpos(fp, &position_before);
    
    while(scheck >= 0) {
        
        if( mGridDumpOn ) {
            scheck = fscanf(fp, "%i %i %e %s %s", &dump, &cycle, &time, &fname, &gname);
            
            if(scheck >= 0) {
                mDumpIDArray.push_back(dump);
                mCycleArray.push_back(cycle);
                mTimeArray.push_back(time);
                mDataFileNameArray.push_back(string(fname));
                mGridFileNameArray.push_back(string(gname));
                
                if(mMPI.getProcID() == 0 && mVerbose) printf("%i %i %e %s %s\n", dump, cycle, time, &fname, &gname);
            }
            
        } else {
            scheck = fscanf(fp, "%i %i %e %s", &dump, &cycle, &time, &fname);
            
            if(scheck >= 0) {
                mDumpIDArray.push_back(dump);
                mCycleArray.push_back(cycle);
                mTimeArray.push_back(time);
                mDataFileNameArray.push_back(string(fname));
                
                if(mMPI.getProcID() == 0 && mVerbose) printf("%i %i %e %s\n", dump, cycle, time, &fname);
            }
        }
    }
}


//====================================================================
// READ AND STORE DATA
//====================================================================
void PPCLASS:: readData(int dumpID)
{
    clearFields();
    
    // define HDF5 variables and arrays
    hid_t    aid1,attr1;
    hid_t    memspaceID,fileID,spaceID,datasetID,datatypeID;
    hid_t    status,filespaceID,cparms;
    hsize_t  count[1],dims[1],dims_chunk[1],dims_nodes[1],dims_vector[1],dims_parts[1];
    hsize_t  offset[1];
    
    //------------------------------------------------------------------
    // specify which data dump cycle to read, and set up input arrays
    //------------------------------------------------------------------
    int rank, rank_chunk;
    char* fieldtoread;
    string sfilename;
    
    int scalarFieldOffset  = 0;
    int vectorFieldOffset  = 0;
    int nodePositionOffset = 0;
    vector<int> particleScalarFieldOffset(mNumParticleTypes,0);
    vector<int> particleVectorFieldOffset(mNumParticleTypes,0);
    
    mScalarFields.resize(mNumScalarFields);
    mVectorFields.resize(mNumVectorFields);
    mParticleScalarFields.resize(mNumParticleTypes);
    mParticleVectorFields.resize(mNumParticleTypes);
    
    // ndir1 and ndir2 specify the lower and upper bound of the
    // processor/directory domain data to be read.
    // The combination ndir1=0, ndir2=mNumDomains specifies to
    // read in the entire data set
    int ndir1 = 0;
    int ndir2 = mNumDomains;
    int nSets = mNumDomains/mMPI.getNumProcs();
    
    if(mNumDomains % mMPI.getNumProcs() != 0) {
        if(mMPI.getProcID() == 0) {
            cout << endl;
            cout << "****WARNING: number of directories is not evenly  " << endl;
            cout << "             divisible by the number of processors" << endl;
            cout << "             change the number of processors" << endl;
        }
        exit(0);
    }
    
    // distribute equally the load to read data across processors
    // each processor reads in nSets directories
    if(mMPI.getNumProcs() > 1) {
        ndir1 = nSets*mMPI.getProcID();
        ndir2 = ndir1 + nSets;
    }
    
    mNumZones = 0;
    for(int idir = ndir1; idir < ndir2; idir++) {
        sfilename = mRootDirName + mDirNameArray[idir]
        + string("/") + mDataFileNameArray[dumpID];
        
        const char *filename = sfilename.data();
        const char *dataname = "Cosmos++";
        if(mMPI.getProcID() == 0 && mVerbose) {
            cout << " " << endl;
            printf("Reading data from file %s\n", filename);
        }
        
        
        //------------------------------------------------------------------
        // open HDF5 dump files, retrieve data and chunk sizes (total zones)
        //------------------------------------------------------------------
        fileID    = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
        datasetID = H5Dopen(fileID, dataname, H5P_DEFAULT);
        
        filespaceID = H5Dget_space(datasetID);
        rank        = H5Sget_simple_extent_ndims(filespaceID);
        status      = H5Sget_simple_extent_dims(filespaceID, dims, NULL);
        if(mMPI.getProcID() == 0 && mVerbose) printf("dataset rank %d, dimension %lu \n",
                                                     rank, (unsigned long)(dims[0]) );
        
        cparms = H5Dget_create_plist(datasetID);
        if(H5D_CHUNKED == H5Pget_layout(cparms)) {
            rank_chunk = H5Pget_chunk(cparms, 1, dims_chunk);
            if(mMPI.getProcID() == 0 && mVerbose) printf("chunk rank %d, dimension %lu\n",
                                                         rank_chunk, (unsigned long)(dims_chunk[0]) );
        }
        
        
        //------------------------------------------------------------------
        // retrieve internal zone size attribute
        //------------------------------------------------------------------
        int numInternalZones;
        attr1  = H5Aopen_name(datasetID, "Number of Internal Zones");
        status = H5Aread(attr1, H5T_NATIVE_INT, &numInternalZones);
        if(mMPI.getProcID() == 0 && mVerbose) printf("Number of Internal Zones = %d\n",
                                                     numInternalZones);
        status = H5Aclose(attr1);
        
        mNumZones += numInternalZones;
        
        
        //------------------------------------------------------------------
        // retrieve particle type size attribute
        //------------------------------------------------------------------
        if(mNumParticleTypes > 0) {
            attr1  = H5Aopen_name(datasetID, "Number of Particle Types");
            status = H5Aread(attr1, H5T_NATIVE_INT, &mNumParticleTypes);
            if(mMPI.getProcID() == 0 && mVerbose) printf("Number of Particle Types = %d\n",
                                                         mNumParticleTypes);
            status = H5Aclose(attr1);
            for(int type = 0; type < mNumParticleTypes; type++) {
                const char *partName = mParticleTypeNames[type].data();
                attr1  = H5Aopen_name(datasetID, partName);
                status = H5Aread(attr1, H5T_NATIVE_INT, &mNumParticles[type]);
                if(mMPI.getProcID() == 0 && mVerbose) printf("Number of Particles of Type %s = %d\n",
                                                             partName, mNumParticles[type]);
                status = H5Aclose(attr1);
            }
        }
        
        //------------------------------------------------------------------
        // Read all field data by chunk (chunk size = number of zones).
        // But only the internal zones are stored:
        //     domain and physical boundaries are ignored.
        //------------------------------------------------------------------
        int ndims = mDimension;
        int nsize = numInternalZones;
        
        dims_nodes[0]    = mNumNodesInZone*ndims*dims_chunk[0];
        dims_vector[0]   = ndims*dims_chunk[0];
        int maxArraySize = dims_nodes[0];
        
        double *chunk_data;
        chunk_data = new double[maxArraySize];
        
        count [0]  = 0;
        offset[0]  = 0;
        
        
        // scalar fields
        // -------------
        int icount = 0;
        for(int ifield = 0; ifield < mNumScalarFields; ifield++) {
            memspaceID = H5Screate_simple(rank_chunk, dims_chunk, NULL);
            offset[0] += count[0];
            count[0]   = dims_chunk[0];
            status = H5Sselect_hyperslab(filespaceID, H5S_SELECT_SET, offset,
                                         NULL, count, NULL);
            status = H5Dread(datasetID, H5T_NATIVE_DOUBLE, memspaceID, filespaceID,
                             H5P_DEFAULT, chunk_data);
            
            mScalarFields[icount].resize( scalarFieldOffset + nsize );
            for(int n = 0; n < nsize; n++) {
                mScalarFields[icount][scalarFieldOffset + n] = chunk_data[n];
            }
            icount++;
        }
        
        
        // vector fields
        // ---------------
        icount = 0;
        for(int ifield = 0; ifield < mNumVectorFields; ifield++) {
            memspaceID = H5Screate_simple(rank_chunk, dims_vector, NULL);
            offset[0]  += count[0];
            count[0]    = dims_vector[0];
            status = H5Sselect_hyperslab(filespaceID, H5S_SELECT_SET, offset,
                                         NULL, count, NULL);
            status = H5Dread(datasetID, H5T_NATIVE_DOUBLE, memspaceID, filespaceID,
                             H5P_DEFAULT, chunk_data);
            
            mVectorFields[icount].resize( vectorFieldOffset + ndims*nsize );
            
            int ncnt = 0;
            for(int n = 0; n < nsize; n++) {
                for(int i = 0; i < mDimension; i++) {
                    mVectorFields[icount][vectorFieldOffset+ncnt] = chunk_data[ncnt];
                    ncnt += 1;
                }
            }
            icount++;
        }
        
        
        // Particle data
        for(int type = 0; type < mNumParticleTypes; type++) {
            mParticleScalarFields[type].resize(mNumParticleScalarFields[type]);
            mParticleVectorFields[type].resize(mNumParticleVectorFields[type]);

            // particle scalar fields
            // -------------
            icount = 0;
            for(int ifield = 0; ifield < mNumParticleScalarFields[type]; ifield++) {
                dims_parts[0] = mNumParticles[type];
                double *part_data;
                part_data = new double[mNumParticles[type]];
                
                memspaceID = H5Screate_simple(rank_chunk, dims_parts, NULL);
                offset[0] += count[0];
                count[0]   = dims_parts[0];
                status = H5Sselect_hyperslab(filespaceID, H5S_SELECT_SET, offset,
                                             NULL, count, NULL);
                status = H5Dread(datasetID, H5T_NATIVE_DOUBLE, memspaceID, filespaceID,
                                 H5P_DEFAULT, part_data);
                
                mParticleScalarFields[type][icount].resize( particleScalarFieldOffset[type] + mNumParticles[type] );
                for(int n = 0; n < mNumParticles[type]; n++) {
                    mParticleScalarFields[type][icount][particleScalarFieldOffset[type] + n] = part_data[n];
                }
                icount++;
                delete[] part_data;
            }
            
            
            // particle vector fields
            // ---------------
            icount = 0;
            for(int ifield = 0; ifield < mNumParticleVectorFields[type]; ifield++) {
                dims_parts[0] = ndims*mNumParticles[type];
                double *part_data;
                part_data = new double[ndims*mNumParticles[type]];
                
                memspaceID = H5Screate_simple(rank_chunk, dims_parts, NULL);
                offset[0]  += count[0];
                count[0]    = dims_parts[0];
                status = H5Sselect_hyperslab(filespaceID, H5S_SELECT_SET, offset,
                                             NULL, count, NULL);
                status = H5Dread(datasetID, H5T_NATIVE_DOUBLE, memspaceID, filespaceID,
                                 H5P_DEFAULT, part_data);
                
                mParticleVectorFields[type][icount].resize( particleVectorFieldOffset[type] + ndims*mNumParticles[type] );
                
                int ncnt = 0;
                for(int n = 0; n < mNumParticles[type]; n++) {
                    for(int i = 0; i < mDimension; i++) {
                        mParticleVectorFields[type][icount][particleVectorFieldOffset[type]+ncnt] = part_data[ncnt];
                        ncnt += 1;
                    }
                }
                icount++;
                delete[] part_data;
            }
            particleScalarFieldOffset[type]  += mNumParticles[type];
            particleVectorFieldOffset[type]  += ndims*mNumParticles[type];
        }
        
        
        if( !mGridDumpOn ) {
            // node positions
            // --------------
            memspaceID = H5Screate_simple(rank_chunk, 4, NULL);
            offset[0]  += count[0];
            count[0]    = dims_nodes[0];
            status = H5Sselect_hyperslab(filespaceID, H5S_SELECT_SET, offset,
                                         NULL, count, NULL);
            status = H5Dread(datasetID, H5T_NATIVE_DOUBLE, memspaceID,
                             filespaceID, H5P_DEFAULT, chunk_data);
            
            int ncnt = 0;
            mNodePosition.resize( nodePositionOffset + ndims*mNumNodesInZone*nsize );
            for(int n = 0; n < nsize; n++) {
                for(int nnodes = 0; nnodes < mNumNodesInZone; nnodes++) {
                    for(int i = 0; i < mDimension; i++) {
                        mNodePosition[nodePositionOffset+ncnt] = chunk_data[ncnt];
                        ncnt += 1;
                    }
                }
            }
            nodePositionOffset += mNumNodesInZone*ndims*nsize;
        }
        
        
        // update shift counters
        // ---------------------
        scalarFieldOffset  += nsize;
        vectorFieldOffset  += ndims*nsize;
        
        
        //------------------------------------------------------------------
        // clean up
        //------------------------------------------------------------------
        delete[] chunk_data;
        status = H5Pclose(cparms);
        status = H5Dclose(datasetID);
        status = H5Sclose(filespaceID);
        status = H5Sclose(memspaceID);
        status = H5Fclose(fileID);
        
        
        if( mGridDumpOn ) {
            //------------------------------------------------------------------
            // Grid data
            //------------------------------------------------------------------
            sfilename = mRootDirName + mDirNameArray[idir]
            + string("/") + mGridFileNameArray[dumpID];
            
            const char *filename2 = sfilename.data();
            if(mMPI.getProcID() == 0 && mVerbose) {
                cout << " " << endl;
                printf("Reading grid data from file %s\n", filename2);
            }
            
            //------------------------------------------------------------------
            // open HDF5 dump files, retrieve data and chunk sizes (total zones)
            //------------------------------------------------------------------
            fileID    = H5Fopen(filename2, H5F_ACC_RDONLY, H5P_DEFAULT);
            datasetID = H5Dopen(fileID, dataname, H5P_DEFAULT);
            
            filespaceID = H5Dget_space(datasetID);
            rank        = H5Sget_simple_extent_ndims(filespaceID);
            status      = H5Sget_simple_extent_dims(filespaceID, dims, NULL);
            if(mMPI.getProcID() == 0 && mVerbose) printf("dataset rank %d, dimension %lu \n",
                                                         rank, (unsigned long)(dims[0]) );
            
            cparms = H5Dget_create_plist(datasetID);
            if(H5D_CHUNKED == H5Pget_layout(cparms)) {
                rank_chunk = H5Pget_chunk(cparms, 1, dims_chunk);
                if(mMPI.getProcID() == 0 && mVerbose) printf("chunk rank %d, dimension %lu\n",
                                                             rank_chunk, (unsigned long)(dims_chunk[0]) );
            }
            
            chunk_data = new double[maxArraySize];
            
            count [0]  = 0;
            offset[0]  = 0;
            
            // node positions
            // --------------
            memspaceID = H5Screate_simple(rank_chunk, dims_nodes, NULL);
            offset[0]  += count[0];
            count[0]    = dims_nodes[0];
            status = H5Sselect_hyperslab(filespaceID, H5S_SELECT_SET, offset,
                                         NULL, count, NULL);
            status = H5Dread(datasetID, H5T_NATIVE_DOUBLE, memspaceID,
                             filespaceID, H5P_DEFAULT, chunk_data);
            
            int ncnt = 0;
            mNodePosition.resize( nodePositionOffset + ndims*mNumNodesInZone*nsize );
            for(int n = 0; n < nsize; n++) {
                for(int nnodes = 0; nnodes < mNumNodesInZone; nnodes++) {
                    for(int i = 0; i < mDimension; i++) {
                        mNodePosition[nodePositionOffset+ncnt] = chunk_data[ncnt];
                        ncnt += 1;
                    }
                }
            }
            
            
            // update shift counters
            // ---------------------
            nodePositionOffset += mNumNodesInZone*ndims*nsize;
            
            //------------------------------------------------------------------
            // clean up
            //------------------------------------------------------------------
            delete[] chunk_data;
            status = H5Pclose(cparms);
            status = H5Dclose(datasetID);
            status = H5Sclose(filespaceID);
            status = H5Sclose(memspaceID);
            status = H5Fclose(fileID);
        }
    }
    
    // contruct zone geometry
    setZoneGeometry();
    setBoxDimensions();
}


//====================================================================
// READ AND STORE DATA
//====================================================================
void PPCLASS::readData(int dumpID, int px, int py, int pz,
                       int opx, int opy, int opz)
{
    clearFields();
    
    // define HDF5 variables and arrays
    hid_t    aid1,attr1;
    hid_t    memspaceID,fileID,spaceID,datasetID,datatypeID;
    hid_t    status,filespaceID,cparms;
    hsize_t  count[1],dims[1],dims_chunk[1],dims_nodes[1],dims_vector[1],dims_parts[1];
    hsize_t  offset[1];
    
    //------------------------------------------------------------------
    // specify which data dump cycle to read, and set up input arrays
    //------------------------------------------------------------------
    int rank, rank_chunk;
    char* fieldtoread;
    string sfilename;
    
    int scalarFieldOffset  = 0;
    int vectorFieldOffset  = 0;
    int nodePositionOffset = 0;
    vector<int> particleScalarFieldOffset(mNumParticleTypes,0);
    vector<int> particleVectorFieldOffset(mNumParticleTypes,0);
    
    mScalarFields.resize(mNumScalarFields);
    mVectorFields.resize(mNumVectorFields);
    mParticleScalarFields.resize(mNumParticleTypes);
    mParticleVectorFields.resize(mNumParticleTypes);
    
    // ndir1 and ndir2 specify the lower and upper bound of the
    // processor/directory domain data to be read.
    // The combination ndir1=0, ndir2=mNumDomains specifies to
    // read in the entire data set
    int ndir1 = 0;
    int ndir2 = mNumDomains;
    int nSets = mNumDomains/mMPI.getNumProcs();
    
    if(mNumDomains % mMPI.getNumProcs() != 0) {
        if(mMPI.getProcID() == 0) {
            cout << endl;
            cout << "****WARNING: number of directories is not evenly  " << endl;
            cout << "             divisible by the number of processors" << endl;
            cout << "             change the number of processors" << endl;
        }
        exit(0);
    }
    
    // distribute equally the load to read data across processors
    // each processor reads in nSets directories
    if(mMPI.getNumProcs() > 1) {
        ndir1 = nSets*mMPI.getProcID();
        ndir2 = ndir1 + nSets;
    }
    
    int pxr = opx/px; int pyr = opy/py; int pzr = opz/pz;
    vector<int> mSets(nSets);
    //int mStartID;
    //for(int kp = 0; kp < pz; kp++) {
    //  for(int jp = 0; jp < py; jp++) {
    //    for(int ip = 0; ip < px; ip++) {
    //	mStartID = kp*pzr*opx*opy + jp*pyr*opx + ip*pxr;
    //    }
    //  }
    //}
    int lProcID = (mMPI.getProcID() % (px*py*pz));
    int mStartID = (lProcID/(px*py))*pzr*opx*opy +
    ((lProcID % (px*py))/px)*pyr*opx +
    (lProcID % px)*pxr;
    //cout << lProcID << " " << (lProcID/(px*py))*pzr*opx*opy << endl;
    //cout << ((lProcID % (px*py))/px)*pyr*opx << " " << (lProcID % px)*pxr << endl;
    //cout << mStartID << endl;
    mStartID += (mMPI.getProcID()/(px*py*pz))*opx*opy*opz;
    //cout << mStartID << endl;
    int n = 0;
    for(int kp = 0; kp < pzr; kp++) {
        for(int jp = 0; jp < pyr; jp++) {
            for(int ip = 0; ip < pxr; ip++) {
                mSets[n] = kp*opx*opy + jp*opx + ip + mStartID;
                n++;
            }
        }
    }
    mNumZones = 0;
    for(int idir = 0; idir < nSets; idir++) {
        sfilename = mRootDirName + mDirNameArray[mSets[idir] ]
        + string("/") + mDataFileNameArray[dumpID];
        
        const char *filename = sfilename.data();
        const char *dataname = "Cosmos++";
        if(mMPI.getProcID() == 0 && mVerbose) {
            cout << " " << endl;
            printf("Reading data from file %s\n", filename);
        }
        
        
        //------------------------------------------------------------------
        // open HDF5 dump files, retrieve data and chunk sizes (total zones)
        //------------------------------------------------------------------
        fileID    = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
        datasetID = H5Dopen(fileID, dataname, H5P_DEFAULT);
        
        filespaceID = H5Dget_space(datasetID);
        rank        = H5Sget_simple_extent_ndims(filespaceID);
        status      = H5Sget_simple_extent_dims(filespaceID, dims, NULL);
        if(mMPI.getProcID() == 0 && mVerbose) printf("dataset rank %d, dimension %lu \n",
                                                     rank, (unsigned long)(dims[0]) );
        
        cparms = H5Dget_create_plist(datasetID);
        if(H5D_CHUNKED == H5Pget_layout(cparms)) {
            rank_chunk = H5Pget_chunk(cparms, 1, dims_chunk);
            if(mMPI.getProcID() == 0 && mVerbose) printf("chunk rank %d, dimension %lu\n",
                                                         rank_chunk, (unsigned long)(dims_chunk[0]) );
        }
        
        
        //------------------------------------------------------------------
        // retrieve internal zone size attribute
        //------------------------------------------------------------------
        int numInternalZones;
        attr1  = H5Aopen_name(datasetID, "Number of Internal Zones");
        status = H5Aread(attr1, H5T_NATIVE_INT, &numInternalZones);
        if(mMPI.getProcID() == 0 && mVerbose) printf("Number of Internal Zones = %d\n",
                                                     numInternalZones);
        status = H5Aclose(attr1);
        
        mNumZones += numInternalZones;
        
        
        //------------------------------------------------------------------
        // retrieve number of particle types and
        // number of particles of each type
        //------------------------------------------------------------------
        if(mNumParticleTypes > 0) {
            attr1  = H5Aopen_name(datasetID, "Number of Particle Types");
            status = H5Aread(attr1, H5T_NATIVE_INT, &mNumParticleTypes);
            if(mMPI.getProcID() == 0 && mVerbose) printf("Number of Particle Types = %d\n",
                                                         mNumParticleTypes);
            status = H5Aclose(attr1);
            for(int type = 0; type < mNumParticleTypes; type++) {
                const char *partName = mParticleTypeNames[type].data();
                attr1  = H5Aopen_name(datasetID, partName);
                status = H5Aread(attr1, H5T_NATIVE_INT, &mNumParticles[type]);
                if(mMPI.getProcID() == 0 && mVerbose) printf("Number of Particles of Type %s = %d\n",
                                                             partName, mNumParticles[type]);
                status = H5Aclose(attr1);
            }
        }
        
        
        //------------------------------------------------------------------
        // Read all field data by chunk (chunk size = number of zones).
        // But only the internal zones are stored:
        //     domain and physical boundaries are ignored.
        //------------------------------------------------------------------
        int ndims = mDimension;
        int nsize = numInternalZones;
        
        dims_nodes[0]    = mNumNodesInZone*ndims*dims_chunk[0];
        dims_vector[0]   = ndims*dims_chunk[0];
        int maxArraySize = dims_nodes[0];
        
        double *chunk_data;
        chunk_data = new double[maxArraySize];
        
        count [0]  = 0;
        offset[0]  = 0;
        
        
        // scalar fields
        // -------------
        int icount = 0;
        for(int ifield = 0; ifield < mNumScalarFields; ifield++) {
            memspaceID = H5Screate_simple(rank_chunk, dims_chunk, NULL);
            offset[0] += count[0];
            count[0]   = dims_chunk[0];
            status = H5Sselect_hyperslab(filespaceID, H5S_SELECT_SET, offset,
                                         NULL, count, NULL);
            status = H5Dread(datasetID, H5T_NATIVE_DOUBLE, memspaceID, filespaceID,
                             H5P_DEFAULT, chunk_data);
            
            mScalarFields[icount].resize( scalarFieldOffset + nsize );
            for(int n = 0; n < nsize; n++) {
                mScalarFields[icount][scalarFieldOffset + n] = chunk_data[n];
            }
            icount++;
        }
        
        
        // vector fields
        // ---------------
        icount = 0;
        for(int ifield = 0; ifield < mNumVectorFields; ifield++) {
            memspaceID = H5Screate_simple(rank_chunk, dims_vector, NULL);
            offset[0]  += count[0];
            count[0]    = dims_vector[0];
            status = H5Sselect_hyperslab(filespaceID, H5S_SELECT_SET, offset,
                                         NULL, count, NULL);
            status = H5Dread(datasetID, H5T_NATIVE_DOUBLE, memspaceID, filespaceID,
                             H5P_DEFAULT, chunk_data);
            
            mVectorFields[icount].resize( vectorFieldOffset + ndims*nsize );
            
            int ncnt = 0;
            for(int n = 0; n < nsize; n++) {
                for(int i = 0; i < mDimension; i++) {
                    mVectorFields[icount][vectorFieldOffset+ncnt] = chunk_data[ncnt];
                    ncnt += 1;
                }
            }
            icount++;
        }
        
        
        // Particle data
        for(int type = 0; type < mNumParticleTypes; type++) {
            mParticleScalarFields[type].resize(mNumParticleScalarFields[type]);
            mParticleVectorFields[type].resize(mNumParticleVectorFields[type]);
            
            // particle scalar fields
            // -------------
            icount = 0;
            for(int ifield = 0; ifield < mNumParticleScalarFields[type]; ifield++) {
                dims_parts[0] = mNumParticles[type];
                double *part_data;
                part_data = new double[mNumParticles[type]];
                
                memspaceID = H5Screate_simple(rank_chunk, dims_parts, NULL);
                offset[0] += count[0];
                count[0]   = dims_parts[0];
                status = H5Sselect_hyperslab(filespaceID, H5S_SELECT_SET, offset,
                                             NULL, count, NULL);
                status = H5Dread(datasetID, H5T_NATIVE_DOUBLE, memspaceID, filespaceID,
                                 H5P_DEFAULT, part_data);
                
                mParticleScalarFields[type][icount].resize( particleScalarFieldOffset[type] + mNumParticles[type] );
                for(int n = 0; n < mNumParticles[type]; n++) {
                    mParticleScalarFields[type][icount][particleScalarFieldOffset[type] + n] = part_data[n];
                }
                icount++;
                delete[] part_data;
            }
            
            
            // particle vector fields
            // ---------------
            icount = 0;
            for(int ifield = 0; ifield < mNumParticleVectorFields[type]; ifield++) {
                dims_parts[0] = ndims*mNumParticles[type];
                double *part_data;
                part_data = new double[ndims*mNumParticles[type]];
                
                memspaceID = H5Screate_simple(rank_chunk, dims_parts, NULL);
                offset[0]  += count[0];
                count[0]    = dims_parts[0];
                status = H5Sselect_hyperslab(filespaceID, H5S_SELECT_SET, offset,
                                             NULL, count, NULL);
                status = H5Dread(datasetID, H5T_NATIVE_DOUBLE, memspaceID, filespaceID,
                                 H5P_DEFAULT, part_data);
                
                mParticleVectorFields[type][icount].resize( particleVectorFieldOffset[type] + ndims*mNumParticles[type] );
                
                int ncnt = 0;
                for(int n = 0; n < mNumParticles[type]; n++) {
                    for(int i = 0; i < mDimension; i++) {
                        mParticleVectorFields[type][icount][particleVectorFieldOffset[type]+ncnt] = part_data[ncnt];
                        ncnt += 1;
                    }
                }
                icount++;
                delete[] part_data;
            }
            particleScalarFieldOffset[type]  += mNumParticles[type];
            particleVectorFieldOffset[type]  += ndims*mNumParticles[type];
        }
        
        
        if( !mGridDumpOn ) {
            // node positions
            // --------------
            memspaceID = H5Screate_simple(rank_chunk, dims_nodes, NULL);
            offset[0]  += count[0];
            count[0]    = dims_nodes[0];
            status = H5Sselect_hyperslab(filespaceID, H5S_SELECT_SET, offset,
                                         NULL, count, NULL);
            status = H5Dread(datasetID, H5T_NATIVE_DOUBLE, memspaceID,
                             filespaceID, H5P_DEFAULT, chunk_data);
            
            int ncnt = 0;
            mNodePosition.resize( nodePositionOffset + ndims*mNumNodesInZone*nsize );
            for(int n = 0; n < nsize; n++) {
                for(int nnodes = 0; nnodes < mNumNodesInZone; nnodes++) {
                    for(int i = 0; i < mDimension; i++) {
                        mNodePosition[nodePositionOffset+ncnt] = chunk_data[ncnt];
                        ncnt += 1;
                    }
                }
            }
            nodePositionOffset += mNumNodesInZone*ndims*nsize;
        }
        
        
        // update shift counters
        // ---------------------
        scalarFieldOffset  += nsize;
        vectorFieldOffset  += ndims*nsize;
        
        
        //------------------------------------------------------------------
        // clean up
        //------------------------------------------------------------------
        delete[] chunk_data;
        status = H5Pclose(cparms);
        status = H5Dclose(datasetID);
        status = H5Sclose(filespaceID);
        status = H5Sclose(memspaceID);
        status = H5Fclose(fileID);
        
        
        if( mGridDumpOn ) {
            sfilename = mRootDirName + mDirNameArray[mSets[idir] ]
            + string("/") + mGridFileNameArray[dumpID];
            
            const char *filename2 = sfilename.data();

            //--------------            if(mMPI.getProcID() == 0 && mVerbose) {
                cout << " " << endl;
                printf("Reading grid data from file %s\n", filename2);
            }
            ----------------------------------------------------
            // open HDF5 dump files, retrieve data and chunk sizes (total zones)
            //------------------------------------------------------------------
            fileID    = H5Fopen(filename2, H5F_ACC_RDONLY, H5P_DEFAULT);
            datasetID = H5Dopen(fileID, dataname, H5P_DEFAULT);
            
            filespaceID = H5Dget_space(datasetID);
            rank        = H5Sget_simple_extent_ndims(filespaceID);
            status      = H5Sget_simple_extent_dims(filespaceID, dims, NULL);
            if(mMPI.getProcID() == 0 && mVerbose) printf("dataset rank %d, dimension %lu \n",
                                                         rank, (unsigned long)(dims[0]) );
            
            cparms = H5Dget_create_plist(datasetID);
            if(H5D_CHUNKED == H5Pget_layout(cparms)) {
                rank_chunk = H5Pget_chunk(cparms, 1, dims_chunk);
                if(mMPI.getProcID() == 0 && mVerbose) printf("chunk rank %d, dimension %lu\n", 
                                                             rank_chunk, (unsigned long)(dims_chunk[0]) );
            }
            
            chunk_data = new double[maxArraySize];
            
            count [0]  = 0;
            offset[0]  = 0;
            
            // node positions 
            // --------------
            memspaceID = H5Screate_simple(rank_chunk, dims_nodes, NULL);
            offset[0]  += count[0];
            count[0]    = dims_nodes[0];
            status = H5Sselect_hyperslab(filespaceID, H5S_SELECT_SET, offset, 
                                         NULL, count, NULL);
            status = H5Dread(datasetID, H5T_NATIVE_DOUBLE, memspaceID, 
                             filespaceID, H5P_DEFAULT, chunk_data);
            
            int ncnt = 0;
            mNodePosition.resize( nodePositionOffset + ndims*mNumNodesInZone*nsize );
            for(int n = 0; n < nsize; n++) {
                for(int nnodes = 0; nnodes < mNumNodesInZone; nnodes++) {
                    for(int i = 0; i < mDimension; i++) {
                        mNodePosition[nodePositionOffset+ncnt] = chunk_data[ncnt];
                        ncnt += 1;
                    }
                }
            }
            
            // update shift counters
            // ---------------------
            nodePositionOffset += mNumNodesInZone*ndims*nsize;
            
            //------------------------------------------------------------------
            // clean up
            //------------------------------------------------------------------
            delete[] chunk_data;
            status = H5Pclose(cparms);
            status = H5Dclose(datasetID);
            status = H5Sclose(filespaceID);
            status = H5Sclose(memspaceID);
            status = H5Fclose(fileID);
        }
    }
    
    // contruct zone geometry
    setZoneGeometry();
    setBoxDimensions();
}


//====================================================================
// RETURN NODE POSITIONS OF SINGLE CELL
//====================================================================
void PPCLASS::getNodePosition(int izone, vector<double> &nodepos)
{
    nodepos.clear();
    int offset = mNumNodesInZone*mDimension*izone;
    
    int ncnt = 0;
    for(int n = 0; n < mNumNodesInZone; n++) {
        for(int i = 0; i < mDimension; i++) {
            nodepos.push_back(mNodePosition[offset+ncnt]);
            ncnt += 1;
        }
    }
}


//====================================================================
// RETURN INTEGER LIST OF NEIGHBOR ZONE IDS
//====================================================================
void PPCLASS::getNeighborList(int izone, vector<int> &neighbors)     
{
    double tol  = 1.0e-8;
    
    double dist;
    double xmin,xmax,ymin,ymax,zmin,zmax;
    
    xmin = mLocalBoxDimensions[0];
    xmax = mLocalBoxDimensions[1];
    
    if(mDimension > 1) {
        ymin = mLocalBoxDimensions[2];
        ymax = mLocalBoxDimensions[3];
    }
    
    if(mDimension > 2) {
        zmin = mLocalBoxDimensions[4];
        zmax = mLocalBoxDimensions[5];
    }
    
    neighbors.clear();
    vector<double> mynodepos;
    getNodePosition(izone, mynodepos);
    
    for(int i = 0; i < mNumZones; i++) {
        if(i != izone) {
            vector<double> nodepos;
            getNodePosition(i, nodepos);
            int n1 = 0;
            bool addedZone = false;
            while( !addedZone && n1 < nodepos.size() ) {
                int n2 = 0;
                bool isClose = false;
                while( !isClose && n2 < mynodepos.size() ) {
                    if(mDimension == 1) {
                        dist = (fabs(nodepos[n1+0] - mynodepos[n2+0]))/(xmax-xmin);
                        
                    } else if(mDimension == 2) {
                        dist = (fabs(nodepos[n1+0] - mynodepos[n2+0]))/(xmax-xmin)
                        + (fabs(nodepos[n1+1] - mynodepos[n2+1]))/(ymax-ymin);
                        
                    } else {
                        dist = (fabs(nodepos[n1+0] - mynodepos[n2+0]))/(xmax-xmin)
                        + (fabs(nodepos[n1+1] - mynodepos[n2+1]))/(ymax-ymin)
                        + (fabs(nodepos[n1+2] - mynodepos[n2+2]))/(zmax-zmin);
                    }
                    
                    if(fabs(dist) < tol) {
                        isClose = true;
                    }
                    n2 += mDimension;
                }
                
                if(isClose) {
                    addedZone = true;
                    neighbors.push_back(i);
                }
                
                n1 += mDimension;
            }
        } // if(izone)
    } //for(int i)
    
    /*
     // check size of neighbor list for consistency
     // (currently valid only for single level grids)
     bool isConsistent = false;
     int size = neighbors.size();
     
     if(mDimension == 1) {
     if(size == 1 || size == 2) isConsistent = true;
     
     } else if(mDimension == 2) {
     if(size == 3 || size == 5 || size == 8) isConsistent = true;
     
     } else {
     if(size ==  7 || size == 11 || 
     size == 17 || size == 26 ) isConsistent = true;
     }
     
     if(!isConsistent) {
     if(mMPI.getProcID() == 0) {
     cout << endl;
     cout << "***WARNING: neighbor list size is not correct" << endl;
     }
     exit(0);
     }
     */
}



//====================================================================
void PPCLASS::writeDataFile(string fileName, char const *append, int n, ...)
{
    va_list ap;
    va_start(ap, n);
    double val;
    
    FILE *OUT;
    OUT = fopen(fileName.data(), append);
    if(OUT == NULL) {
        printf("problem opening output file in PostProcess\n");  
        exit(0); 
    }
    
    if(n>0) {
        val = va_arg(ap, double);
        fprintf(OUT,"%12.5e", val);
    }
    for(n-=1;n>0;n--) {
        val = va_arg(ap, double);
        fprintf(OUT,"\t%12.5e", val);
    }
    fprintf(OUT,"\n");
    fclose(OUT);
}

//==================================================================
void PPCLASS::writeDataFileArray(string fileName, char const *append, int maxZone, double time, vector<double> data)
{
    double val;
    FILE *OUT;
    OUT = fopen(fileName.data(), append);
    if(OUT == NULL) {
        printf("problem opening output file in PostProcess\n");
        exit(0);
    }
    int i = 0;
    fprintf(OUT, "%12.5e", time);
    for(i = 0; i < maxZone; i++) {
        fprintf(OUT,"\t%12.5e", data[i]);
    }
    fprintf(OUT,"\n");
    fclose(OUT);
}

