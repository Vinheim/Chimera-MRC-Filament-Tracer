/**
 *	(Date: 02/20/2019)
 *	Declaration of a class to support interaction with MRC files containing three-dimensional grids of voxels, 
 *	each of which has a value corresponding to electron density or electric potential. The functions of this class will
 *	allow the user to properly graph and read data regarding actin filaments of hair cell stereocilia.
 *	@author John Hessefort
 */

#ifndef __MRC_H__
#define __MRC_H__

#include <cstdio>
#include <vector>
#include <string>

using namespace std;

struct Coordinate
{
	float x;
	float y;
	float z;
	
	/**
	 *	Default Constructor
	 *	Initialize all coordinate values to be zero.
	 */
	Coordinate() { x = y = z = 0; }

	/**
	 *	Multiple Argument Constructor
	 *	@param C1 : The value used to initialize this coordinate's X-Plane value.
	 *	@param C2 : The value used to initialize this coordinate's Y-Plane value.
	 *	@param C3 : The value used to initialize this coordinate's Z-Plane value.
	 */
	Coordinate(float C1, float C2, float C3)
		: x(C1), y(C2), z(C3) {}
};

struct Index
{
	int x;
	int y;
	int z;

	/**
	 *	Default Constructor
	 *	Initialize all index values to be zero.
	 */
	Index() { x = y = z = 0; }

	/**
	 *	Multiple Argument Constructor
	 *	@param I1: The value used to initialize this index's first index value.
	 *	@param I2: The value used to initialize this index's second index value.
	 *	@param I3: The value used to initialize this index's third index value.
	 */
	Index(int I1, int I2, int I3)
		: x(I1), y(I2), z(I3) {}

	void operator=(const Index& ID)
	{
		if(this != &ID)
		{
			x = ID.x;
			y = ID.y;
			z = ID.z;
		}
	}
};

/**
 *	The Color struct contains information related to a color as expressed through the additive RGB color model.
 *	Three colors: red, green, and blue can be added together in various ways to produce a broad assortment of colors.
 */
struct Color
{
	float r;
	float g;
	float b;

	/**
	 *	Default Constructor
	 *	Initialize all color values: red, green, and blue to be zero.
	 */
	Color() { r = g = b = 0; }

	/**
	 *	Multiple Argument Constructor
	 *	@param C1 : The value used to initialize the red member variable for this color.
	 *	@param C2 : The value used to initialize the green member variable for this color.
	 *	@param C3 : The value used to initialize the blue member variable for thie color.
	 *	@precondition: 0 <= C1, C2, C3 <= 255 
	 */
	Color(float C1, float C2, float C3)
		: r(C1), g(C2), b(C3) {}
};

struct FilamentPlane
{
	int first;
	int last;
	int leastEndDistanceFilament;
};

struct FilamentError
{
	int filamentID;
	float error;
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
	void printDensityAtIndex(const Index& index);
}

#endif // Defined __MRC_H__
