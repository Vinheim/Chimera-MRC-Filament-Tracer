#ifndef ____MRC__
#define ____MRC__

#include <stdio.h>
#include <vector>
#include <string>


using namespace std;

struct Coordinate{
    float xCor;
    float yCor;
    float zCor;
    
    Coordinate(){
    }
    
    
    Coordinate(float x, float y, float z){
        xCor = x;
        yCor = y;
        zCor = z;
    }
};

struct Index{
    int xIndex;
    int yIndex;
    int zIndex;
    
    Index(){
    }
    
    Index(int x, int y, int z){
        xIndex = x;
        yIndex = y;
        zIndex = z;
    }
    
};

struct RGBColor{
    
    float r;
    float g;
    float b;
    
    RGBColor(){
    }
    
    RGBColor(int red, int green, int blue){
        r = (float)red / 255.0;
        g = (float)green / 255.0;
        b = (float)blue / 255.0;
    }
    
};

struct FilamentPlane {
    int firstFilament;
    int lastFilament;
    int leastEndDistanceFilament;
};

struct FilamentError{
    int filamentID;
    float error;
};

struct ManualMarker {
    float distanceToClosestSeedPoint;
    int closestSeedPoint;
    
    ManualMarker() {
        distanceToClosestSeedPoint = 100000000;
        closestSeedPoint = -1;
    }
};

enum FilamentDirection {
    Ascending = 1,
    Descending = -1
};

enum FilamentComparisonOrder {
    FilamentComparisonOrderAscending = 1,
    FilamentComparisonOrderDescending = -1
};

class MRC
{

private:
    int mode;//data type
    //0 image: signed 8-bit bytes range -128 to 127
    //1 image: 16-bit halfwords
    //2 image: 32-bit reals
    //3 transform: complex 16-bit integers
    //4 transform: complex 32-bit reals
    //6 image: unsigned 16-bit range 0 to 65535
    int nxStart;//number of first column in map (Defalut = 0)
    int nyStart;//number of first row in map
    int nzStart;//number of first section in map
    int mx;//number of intervals along X
    int my;//number of intervals along Y
    int mz;//number of intervals along Z
    float xLength;//cell dimensions in X (angstrom)
    float yLength;//cell dimensions in Y (angstrom)
    float zLength;//cell dimensions in Z (angstrom)
    float alpha;//cell angle between Y and Z
    float beta;//cell angle between X and Z
    float gamma;//cell angle between X and Y
    int mapc;//axis corresp to cols(1,2,3 for X,Y,Z)
    int mapr;//axis corresp to rows(1,2,3 for X,Y,Z)
    int maps;//axis corresp to sections(1,2,3 for X,Y,Z)
    float dMin;//minimum density value
    float dMax;//maximum density value
    float dMean;//mean density value
    int ispg;//space group number 0 or 1 (default = 0)
    int nsymbt;//number of bytes used for symmetry data (0 or 80)
    int extra[25];//extra space used for anything - 0 by default
    
    char map[4];//character string 'MPA' to identify file type
    int machineStamp;//machine stamp
    float rms;//rms deviation of map from mean density
    int nlabl;//number of labels being used
    char label[10][80];//10 80-character text labels
    //?? make it private later float ***cube;//matrix points of density map
    
    
public:
    
    int nx;//number of columns (fastest changing in map)
    int ny;//number of rows
    int nz;//number of sections (slowest changing in map)
    
    float ***cube;
    float ***pathDensity;
    int ***path;
    
    float xOrigin;//origin in X used for transforms
    float yOrigin;//origin in Y used for transforms
    float zOrigin;//origin in Z used for transforms
    
    double getApixX() const {return xLength/mx;}
    double getApixY() const {return yLength/my;}
    double getApixZ() const {return zLength/mz;}
    
    MRC()//consturctor
    {
    }
    
    MRC(MRC &right);//copy constructor
    
    
    void readMRCFile(string inputMRCFile);//read density map by EMan formate
    void writeDataInMRCFormat(string outputMRCFile);
    void showMRCFileDescription();
    float minDensity();
    float maxDensity();
    float meanDensity();
    
    //Gaussian Matrix
    void initializeGaussianPeakCube();
    void copyGaussianPeakCubeToMRCCube(int yIndex, int range);
    void calculatePeakValue(float &peakValue, int i, int j, float pMX, float pMY, float pSX, float pSZ);
    void findSeedPointsFromMRC();
    void generateSevenPeaksGaussianMatrix(vector < vector <vector<float> > > &allPeakDistribution);
    void findTheMaximumDensityPointAtYCrossSection(int yCrossSection);
    void writeInPDBFormat(string outputFile, vector<Coordinate> coordinates);
    void findCenterOfMass(int yCrossSection, int &xCenter, int &zCenter);
    
    //Tracing
    void readSeedPointsNTraceFilaments(string seedPointsFilename);
    void findClosestManualModels(vector<Coordinate> seedPoints);
    void traceFilamentsWithStepSize(vector<Coordinate> seedpoints, int  maximumXShiftInBoundngBox);
    void traceFilament(Index startVoxel, Index endVoxel, vector <Coordinate> &
                       denseVoxels, int maximumXShiftInBoundngBox, int filamentID);
    
    void getIndexFromCoordinate(float xCoord, int yCoord, float zCoord, Index &seedIndex);
    void printFilamentPath(Index start, Index end , int filamentId);
    void writeDataInCMMFormat(string outputFile, vector<Coordinate> coordinates, string filamentName,int filamentId);
    
    void getCoordinateFromIndex(int xIndex, int yIndex, int zIndex, Coordinate &coordinate);
    int getIndexAtAxisFromCoordinate(float CoordAtAxis, float originAtAxis, float Api,float startAtAxis);
    float getCoordinateAtAxisFromIndex(float IndexAtAxis, float originAtAxis, float Api,float startAtAxis);
    
    void findPath(Index voxel, int stepSize, int yDirection, int &xIndex);
    
    //void findMaxDenseVoxelFromRow(Index &voxel, int xWindowSize);
    void findNextReferenceVoxel(Index refvoxel, Index &newRefVoxel, int stepSize, int xShift, float &avgDensityinBoundingBox);
    void findMaxDenseVoxelFromRow(Index &voxel, int xWindowSizeLeft, int xWindowSizeRight);
    void findMaxDenseVoxelFromTopRow(Index &voxel, int xWindowSizeLeft, int xWindowSizeRight, int boundinBoxSize, float &avgDensityinBoundingBox);
    
    void resetPathDensity(Index voxel, int windowSizeInX);
    void resetDensityOfBoundingBox(Index voxel, int windowSizeInX);
    
    //Smoothing
    void smoothAllFilaments();
    void smoothFilamentUsingInterpolation(vector <Coordinate> &filament);
    void smoothFilamentWithID(int filamentID);
    void smoothFilamentWithFilename(string filamentInputFile, vector<Coordinate> &smoothedCoordinates);
    
    
    void smoothFilamentWithID(int filamentID, int numberOfIteration);
    void smoothFilamentUsingInterpolation(vector <Coordinate> &filament, int numOfIteration);
    
    //others
    void copyFilamentCoordinates(vector<Coordinate> &toFilament, vector <Coordinate>fromFilament);
    void findXCrspndToY(float filament1CoordY, vector<Coordinate> filament, float &filament2CoordX);
    void findZCrspndToY(float filamentCoordY, vector<Coordinate> filament, float &filamentCoordZ);
    bool isPointOnFilament(float filamentCoordY, vector<Coordinate> filament);
    float calculateXcoordDistanceBetweenFilamentEndpoints(int filamentID1, int filamentID2);
    float getNumberOfOverlappingSegments();
    
    //Accuracy
    void readAndCalculateAccuracy();
  //  Retracing
    bool doesOverlap(float xCoord1, float zCoord1, float xCoord2, float zCoord2);
    bool doesDiverge(float xCoord1, float zCoord1, float xCoord2, float zCoord2);
    void findOverlapingFilaments(vector <vector <Coordinate> > filaments);
    void checkDiversion(vector<Coordinate> filamentToCheck, vector <Coordinate> filamentToCompare, bool &doesDivergeNeighbour, float &diversionCoordY);
    void checkOverlapping(vector<Coordinate> filamentToCheck, vector<Coordinate> filamentToCompare, bool &doesOverlapWithNeighbour,float &intersectionCoordY);
    //void findLastCoordinateFromFilamentFile(string filename, Coordinate &lastCoordinate);
    
    
    void findNFixOverlapingFilaments();
    void fixFilament(float startYCoord, vector<Coordinate> filamentToFix, vector<Coordinate> neighbourFilament, int numberOfiteration , int xDirection, int filamentID);
    void fixAllFilament(float startYCoord, vector<Coordinate> filamentToFix, vector<Coordinate> neighbourFilament, int numOfIteration, int xDirection, FilamentComparisonOrder filamentOrder,int filamentID);
    
    void getPathBeforeIntersection(float newStartingY, vector<Coordinate> filament, vector<Coordinate> &pathBeforeIntersection);
    void mergePath(vector<Coordinate> beforeIntersection, vector<Coordinate> afterIntersection, vector<Coordinate> &modifiedPath);
    void goBackwardUntilNoIntersectionOrDiversion(float intersectionCoordY , vector<Coordinate> filamentToCheck, vector<Coordinate> filamentToCompare, float &beginYCoord);
    
    void retraceFilament(Index voxel, Index startVoxel, Index endVoxel,  vector <Coordinate> &denseVoxels, int offsetXLeft, int offsetXRight, int numOfIteration);
    
    void retracingStep2();
    
    
    //Other
    void averageDensityAlongZ(Index index, float &avgDensity);
    void markFirstFilamentOfPlane();

    
    //Read/write
    void writeDataInCMMFormatWithColor(string outputFile, vector<Coordinate> coordinates, string filamentName,int filamentId, RGBColor rgbColor);
    void writeIntersectionDataInCMMFormat(string outputFile, vector<Coordinate> coordinates, string filamentName);
    void readNSortFilamentFromCMM(string filename, vector<Coordinate> &filament);
    void readCoordinateFromCMMFile(vector<Coordinate> &seeds,string inputFilename);
    
	// Get voxel size
	float getCellA0() const; // Returns xLength portion of cell dimensions in angstroms
	float getCellA1() const; // Returns yLength portion of cell dimensions in angstroms
	float getCellA2() const; // Returns zLength portion of cell dimensions in angstroms
	void printVoxelSize() const; // Print total cell dimensions, measured in angstroms, along each axis 

    ~MRC();
 
    void normalizeDensity();
    void setInitialDensity(Index voxel, int windowSizeInX);
    
};

#endif /* defined(____MRC__) */