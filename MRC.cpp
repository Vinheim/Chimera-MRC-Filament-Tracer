#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "MRC.h"

#if defined(_WIN32)
#include <direct.h>
#else 
#include <sys/stat.h>
#include <sys/types.h>
#endif

using namespace std;

#define DENSITY_INFINITY 1000000

#define VOXEL_PATH_DENSITY -100000
#define VOXEL_PARENT -10
#define SOURCE_VOXEL_PATH_DENSITY 0
#define SOURCE_VOXEL_PARENT -1

#pragma mark- MRC
void MRC::readMRCFile(string inputMRCFile)
{
    cout<<"\nReading MRC file: \n " <<inputMRCFile << endl;
    
    fstream inputFile(inputMRCFile.c_str(), ios::in|ios::binary);
    if(!inputFile)
    {
        cout<<"Can not open binary input file "<<inputMRCFile<<endl;
        exit(0);
    }
    
    inputFile.read(reinterpret_cast<char *>(&nx), sizeof(nx));
    inputFile.read(reinterpret_cast<char *>(&ny), sizeof(ny));
    inputFile.read(reinterpret_cast<char *>(&nz), sizeof(nz));
    inputFile.read(reinterpret_cast<char *>(&mode), sizeof(mode));
    inputFile.read(reinterpret_cast<char *>(&nxStart), sizeof(nxStart));
    inputFile.read(reinterpret_cast<char *>(&nyStart), sizeof(nyStart));
    inputFile.read(reinterpret_cast<char *>(&nzStart), sizeof(nzStart));
    inputFile.read(reinterpret_cast<char *>(&mx), sizeof(mx));
    inputFile.read(reinterpret_cast<char *>(&my), sizeof(my));
    inputFile.read(reinterpret_cast<char *>(&mz), sizeof(mz));
    inputFile.read(reinterpret_cast<char *>(&xLength), sizeof(xLength));
    inputFile.read(reinterpret_cast<char *>(&yLength), sizeof(yLength));
    inputFile.read(reinterpret_cast<char *>(&zLength), sizeof(zLength));
    inputFile.read(reinterpret_cast<char *>(&alpha), sizeof(alpha));
    inputFile.read(reinterpret_cast<char *>(&beta), sizeof(beta));
    inputFile.read(reinterpret_cast<char *>(&gamma), sizeof(gamma));
    inputFile.read(reinterpret_cast<char *>(&mapc), sizeof(mapc));
    inputFile.read(reinterpret_cast<char *>(&mapr), sizeof(mapr));
    inputFile.read(reinterpret_cast<char *>(&maps), sizeof(maps));
    inputFile.read(reinterpret_cast<char *>(&dMin), sizeof(dMin));
    inputFile.read(reinterpret_cast<char *>(&dMax), sizeof(dMax));
    inputFile.read(reinterpret_cast<char *>(&dMean), sizeof(dMean));
    inputFile.read(reinterpret_cast<char *>(&ispg), sizeof(ispg));
    inputFile.read(reinterpret_cast<char *>(&nsymbt), sizeof(nsymbt));
    inputFile.read(reinterpret_cast<char *>(extra), 25*sizeof(int));
    inputFile.read(reinterpret_cast<char *>(&xOrigin), sizeof(xOrigin));
    inputFile.read(reinterpret_cast<char *>(&yOrigin), sizeof(yOrigin));
    inputFile.read(reinterpret_cast<char *>(&zOrigin), sizeof(zOrigin));
    inputFile.read(reinterpret_cast<char *>(map), 4*sizeof(char));
    inputFile.read(reinterpret_cast<char *>(&machineStamp), sizeof(machineStamp));
    inputFile.read(reinterpret_cast<char *>(&rms), sizeof(rms));
    inputFile.read(reinterpret_cast<char *>(&nlabl), sizeof(nlabl));
    inputFile.read(reinterpret_cast<char *>(label), 10*80*sizeof(char));
    
    cube = new float **[nx];
   
    for(int i = 0; i < nx; i++)
    {
        cube[i] = new float *[ny];
        for(int j = 0; j < ny; j++)
            cube[i][j] = new float [nz];
    }
    
    pathDensity = new float **[nx];
    path = new int **[nx];
   
    for(int i = 0; i < nx; i++)
    {
        pathDensity[i] = new float *[ny];
        for(int j = 0; j < ny; j++)
            pathDensity[i][j] = new float [nz];
        
        path[i] = new int *[ny];
        for(int j = 0; j < ny; j++)
            path[i][j] = new int [nz];
    }

    for(int i = 0; i < nz; i++) {
        for(int j = 0; j < ny; j++) {
            for(int k = 0; k < nx; k++) {
                inputFile.read(reinterpret_cast<char *>(&cube[k][j][i]), sizeof(float));
                
            }
        }
    }
    
    inputFile.close();
        
    MRC::showMRCFileDescription();
}


/*
 * Write data in MRC format.
 * file is the output file path.
 * set value of all the header field(nx,ny.......,label) before writting to outputMRCFile in binary Format.
 * e.g nx(number of column),ny (num of row),nz(number of section).
 * to calculate dMin, dMax, dMean you may use provided function.
 */

void MRC::writeDataInMRCFormat(string outputMRCFile)
{
    printf("\nWritting MRC data to binary file: %s\n", outputMRCFile.c_str());
  
    fstream outputFile(outputMRCFile.c_str(), ios::out|ios::binary);
    if(!outputFile)
    {
        cout<<"Can not open binary output file "<<outputMRCFile<<endl;
        exit(0);
    }
    
    dMin = minDensity();
    dMax = maxDensity();
    dMean = meanDensity();
    
    outputFile.write(reinterpret_cast<char *>(&nx), sizeof(nx));
    outputFile.write(reinterpret_cast<char *>(&ny), sizeof(ny));
    outputFile.write(reinterpret_cast<char *>(&nz), sizeof(nz));
    outputFile.write(reinterpret_cast<char *>(&mode), sizeof(mode));
    outputFile.write(reinterpret_cast<char *>(&nxStart), sizeof(nxStart));
    outputFile.write(reinterpret_cast<char *>(&nyStart), sizeof(nyStart));
    outputFile.write(reinterpret_cast<char *>(&nzStart), sizeof(nzStart));
    outputFile.write(reinterpret_cast<char *>(&mx), sizeof(mx));
    outputFile.write(reinterpret_cast<char *>(&my), sizeof(my));
    outputFile.write(reinterpret_cast<char *>(&mz), sizeof(mz));
    outputFile.write(reinterpret_cast<char *>(&xLength), sizeof(xLength));
    outputFile.write(reinterpret_cast<char *>(&yLength), sizeof(yLength));
    outputFile.write(reinterpret_cast<char *>(&zLength), sizeof(zLength));
    outputFile.write(reinterpret_cast<char *>(&alpha), sizeof(alpha));
    outputFile.write(reinterpret_cast<char *>(&beta), sizeof(beta));
    outputFile.write(reinterpret_cast<char *>(&gamma), sizeof(gamma));
    outputFile.write(reinterpret_cast<char *>(&mapc), sizeof(mapc));
    outputFile.write(reinterpret_cast<char *>(&mapr), sizeof(mapr));
    outputFile.write(reinterpret_cast<char *>(&maps), sizeof(maps));
    outputFile.write(reinterpret_cast<char *>(&dMin), sizeof(dMin));
    outputFile.write(reinterpret_cast<char *>(&dMax), sizeof(dMax));
    outputFile.write(reinterpret_cast<char *>(&dMean), sizeof(dMean));
    outputFile.write(reinterpret_cast<char *>(&ispg), sizeof(ispg));
    outputFile.write(reinterpret_cast<char *>(&nsymbt), sizeof(nsymbt));
    outputFile.write(reinterpret_cast<char *>(extra), 25*sizeof(int));
    outputFile.write(reinterpret_cast<char *>(&xOrigin), sizeof(xOrigin));
    outputFile.write(reinterpret_cast<char *>(&yOrigin), sizeof(yOrigin));
    outputFile.write(reinterpret_cast<char *>(&zOrigin), sizeof(zOrigin));
    outputFile.write(reinterpret_cast<char *>(map), 4*sizeof(char));
    outputFile.write(reinterpret_cast<char *>(&machineStamp), sizeof(machineStamp));
    outputFile.write(reinterpret_cast<char *>(&rms), sizeof(rms));
    outputFile.write(reinterpret_cast<char *>(&nlabl), sizeof(nlabl));
    outputFile.write(reinterpret_cast<char *>(label), 10*80*sizeof(char));
    
    for(int i = 0; i < nz; i++)
        for(int j = 0; j < ny; j++)
            for(int k = 0; k < nx; k++)
                outputFile.write(reinterpret_cast<char *>(&cube[k][j][i]), sizeof(float));
    outputFile.close();
    
    printf("\nwriting complete\n");
}


/*
 * Print the value of all MRC file header field
 */

void MRC::showMRCFileDescription()
{
    cout<<"MRC File Description: \n";
    cout<<"-------------------------------------------"<<endl;
    cout<<"nx, ny, nz: "<<nx<<","<<ny<<","<<nz<<endl;
    cout<<"mode: "<<mode<<endl;
    cout<<"nxStart, nyStart, nzStart: "<<nxStart<<","<<nyStart<<","<<nzStart<<endl;
    cout<<"mx, my, mz: "<<mx<<","<<my<<","<<mz<<endl;
    cout<<"xLength, yLength, zLength: "<<xLength<<","<<yLength<<","<<zLength<<endl;
    cout<<"alpha, beta, gamma: "<<alpha<<","<<beta<<","<<gamma<<endl;
    cout<<"mapc, mapr, maps: "<<mapc<<","<<mapr<<","<<maps<<endl;
    cout<<"dMin, dMax, dMean: "<<dMin<<","<<dMax<<","<<dMean<<endl;
    cout<<"ispg: "<<ispg<<endl;
    cout<<"nsymbt: "<<nsymbt<<endl;
    for(int i = 0; i < 25; i ++)
        cout<<extra[i]<<" ";
    cout<<endl;
    cout<<"xOrigin, yOrigin, zOrigin: "<<xOrigin<<","<<yOrigin<<","<<zOrigin<<endl;
    cout<<"map: "<<map<<endl;
    cout<<"machineStamp: "<<machineStamp<<endl;
    cout<<"rms: "<<rms<<endl;
    cout<<"nlabl: "<<nlabl<<endl;
    cout<<"-------------------------------------------"<<endl;
    
    cout<<endl;
}

float MRC:: minDensity()
{
    float minimumDensity = DENSITY_INFINITY;
    
    for(int i = 0; i < nz; i++) {
        for(int j = 0; j < ny; j++) {
            for(int k = 0; k < nx; k++) {
                if (minimumDensity > cube[k][j][i]) {
                    minimumDensity = cube[k][j][i];
                }
            }
        }
    }
    
   // printf("Min Density:%f\n", minimumDensity);
    return minimumDensity;
}

float MRC:: maxDensity()
{
    float maximumDensity = -DENSITY_INFINITY;
    
    for(int i = 0; i < nz; i++) {
        for(int j = 0; j < ny; j++) {
            for(int k = 0; k < nx; k++) {
                if (maximumDensity < cube[k][j][i]) {
                    maximumDensity = cube[k][j][i];
                }
            }
        }
    }
    
    cout<<"Max Density: " << maximumDensity << endl;
    return maximumDensity;
}

float MRC:: meanDensity()
{
    float totalDensity = 0;
    
    for(int i = 0; i < nz; i++) {
        for(int j = 0; j < ny; j++) {
            for(int k = 0; k < nx; k++) {
                totalDensity += cube[k][j][i];
            }
        }
    }
    
    float mDensity = (float) totalDensity / (nz * ny * nx);
    cout<<"Mean Density: " <<mDensity << endl;
    return mDensity;
}

float calculateDistance (Coordinate coord1, Coordinate coord2) {
    
    float xDistance = (coord1.xCor - coord2.xCor);
    xDistance *= xDistance;
    
    float yDistance = (coord1.yCor - coord2.yCor);
    yDistance *= yDistance;

    float zDistance = (coord1.zCor - coord2.zCor);
    zDistance *= zDistance;
    
    return  sqrtf(xDistance + yDistance + zDistance);
}

#pragma mark- Index/Coordinate
/**
 * getIndexFromCoordinate() to convert a physical coordinate's x, y, z values into corresponding x, y, and z indices for the cube 3D array/triple pointer variable.
 * Once we have the converted indices x, y, and z, we can access the densities at a certain point through usage of cube[x][y][z].
 */
void MRC::getIndexFromCoordinate(float xCoord, int yCoord, float zCoord, Index &seedIndex)
{
    seedIndex.xIndex = ((xCoord - xOrigin)/ getApixX()) - nxStart;
    seedIndex.yIndex = ((yCoord - yOrigin)/ getApixY()) - nyStart;
    seedIndex.zIndex = ((zCoord - zOrigin)/ getApixZ()) - nzStart;
}

/**
 * getCoordinateFromIndex() to convert the x, y, z indiced of the cube variable to the physical Chimera coordinate's x, y, and z values.
 */
void MRC::getCoordinateFromIndex(int xIndex, int yIndex, int zIndex, Coordinate &coordinate)
{
    coordinate.xCor = (xOrigin +  ((nxStart + xIndex)  * getApixX()));
    coordinate.yCor = (yOrigin +  ((nyStart + yIndex)  * getApixY()));
    coordinate.zCor = (zOrigin +  ((nzStart + zIndex)  * getApixZ()));
}

#pragma mark- Normalize/SetValue
void MRC::normalizeDensity()
{
    float min = minDensity();
   // printf("minimum density: %f\n", min);
    
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                if (min < 0) {// && cube[i][j][k] < 0) {
                    cube[i][j][k] -= min;
                }
                pathDensity[i][j][k] = VOXEL_PATH_DENSITY;
                path[i][j][k] = VOXEL_PARENT;
            }
        }
    }
}

void MRC::setInitialDensity(Index voxel, int windowSizeInX)
{
    resetPathDensity(voxel, windowSizeInX);
    
    pathDensity[voxel.xIndex][voxel.yIndex][voxel.zIndex] = SOURCE_VOXEL_PATH_DENSITY;
    path[voxel.xIndex][voxel.yIndex][voxel.zIndex] = SOURCE_VOXEL_PARENT;
}

void MRC::resetPathDensity(Index voxel, int windowSizeInX)
{
    for(int i = voxel.xIndex - windowSizeInX; i <= voxel.xIndex + windowSizeInX; i++) {
        for(int j = 0; j < ny; j++) {
            for(int k = voxel.zIndex - 1; k <= voxel.zIndex + 1; k++) {
                if ( i < 0 || i >= nx) {
                    continue;
                }
                pathDensity[i][j][k] = VOXEL_PATH_DENSITY;
                path[i][j][k] = VOXEL_PARENT;
            }
        }
    }
}

#pragma mark- Others
void MRC::copyFilamentCoordinates(vector<Coordinate> &toFilament, vector <Coordinate>fromFilament)
{
    toFilament.clear();
    for (int i = 0; i < fromFilament.size(); i++) {
        toFilament.push_back(fromFilament.at(i));
    }
}

#pragma mark -
void MRC::writeDataInCMMFormatWithColor(string outputFile, vector<Coordinate> coordinates, string filamentName,int filamentId, RGBColor rgbColor)
{
    fstream output(outputFile.c_str(), ios::out);
    if(!output)
    {
        cout<<"Can not open output file "<<outputFile<<endl;
     //   cerr << "Error: " << strerror(errno);
        exit(0);
    }
    srand(time(NULL));
    float r = rgbColor.r;
    float g = rgbColor.g;
    float b = rgbColor.b;
    
    
    if (outputFile.find("backtracked") != std::string::npos) {
        
        r = 1.0;//0.7;
        g = 0.2;//0.3;// 1;// 0.3;
        b = 0.2;//1.0;
        
        
        //Temp
        r = 0.0;//0.7;
        g = 1.0;//0.3;// 1;// 0.3;
        b = 0.7;//1.0;
        
        
    }
    else if (outputFile.find("smoothed") != std::string::npos) {
        r = 0.3;//0.7;
        g = 1.0;//1.0;//0.3;// 1;// 0.3;
        b = 0.0;//1.0;
        
        //Temp
        //   r = 0.5;//0.7;
        //  g = 0.3;//1.0;//0.3;// 1;// 0.3;
        // b = 0.1;//1.0;
    }
    else {
        r = 0.6;//0.7;
        g = 0.5;//0.3;// 1;// 0.3;
        b = 1.0;//1.0;
    }
    

    
    
    string  filamentHeader = filamentName + "_" + to_string(filamentId);
    
    output <<"<marker_set name=\"Filament_"<<filamentHeader<<"\">\n";
    
    for (int i = 0; i < coordinates.size(); i++) {
        
        Coordinate coordinate = coordinates.at(i);
        
        output <<"<marker id=\""<<i<<"\""<<" x=\""<<coordinate.xCor<<"\" "<<"y=\""<<coordinate.yCor<<"\" "<<"z=\""<<coordinate.zCor<<"\" "<<"r=\""<<r<<"\" g=\""<<g<<"\" b=\""<<b<<"\" radius=\"17\"/>";
        output<<"\n";
    }
    
    for (int i = 0; i < coordinates.size() - 1; i++) {
        output <<"<link id1=\""<<i<<"\""<<" id2=\""<<i + 1<<"\""<<" r=\""<<r<<"\" g=\""<<g<<"\" b=\""<<b<<"\" radius=\"14\"/>";
        output<<"\n";
    }
    
    output <<"</marker_set>";
    // printf("Writting data in CMM format complete\n");
}

void MRC::readCoordinateFromCMMFile(vector<Coordinate> &seeds,string inputFilename)
{
    fstream inputFile(inputFilename.c_str(), ios::in);
    if(!inputFile)
    {
        cout<<"Can not open binary input file "<<inputFilename<<endl;
      //  exit(1);
    }
    
    string line;
    while (getline(inputFile, line )) {
        //     cout<<line<<"\n";
        string buffer;
        stringstream sStream(line);
        // std::smatch match;
        
        vector <std::string> splitedStrings;
        while (sStream >> buffer) {
            splitedStrings.push_back(buffer);
            // cout<<buffer<<"\n";
        }
        
        Coordinate coordinate;
        for (int i = 0; i < splitedStrings.size(); i++){
            string str = splitedStrings.at(i);
            
            std::size_t foundX = str.find("x=");
            std::size_t foundY = str.find("y=");
            std::size_t foundZ = str.find("z=");
            
            if (foundX!=std::string::npos){
                string temp = str.substr(foundX + 3, str.size() - foundX - 3);
                //  std::cout << "first 'needle' found at: " << foundX << '\n';
                // cout<<temp;
                coordinate.xCor = stof(temp);
                // cout<<"XValue: "<<coordinate.xCor<<"\n";
            }
            
            if (foundY!=std::string::npos){
                string temp = str.substr(foundY + 3, str.size() - foundY - 3);
                // std::cout << "first 'needle' found at: " << foundY << '\n';
                coordinate.yCor = stof(temp);
                //cout<<"YValue: "<<coordinate.yCor<<"\n";
            }
            
            if (foundZ!=std::string::npos){
                string temp = str.substr(foundZ + 3, str.size() - foundZ - 3);
                //  std::cout << "first 'needle' found at: " << foundZ << '\n';
                coordinate.zCor = stof(temp);
                // cout<<"zValue: "<<coordinate.zCor<<"\n";
                seeds.push_back(coordinate);
            }
        }
    }
    
    // printf("Done with the Seeds Reading\n");
}

#pragma mark -
void MRC::writeDataInCMMFormat(string outputFile, vector<Coordinate> coordinates, string filamentName,int filamentId)
{
    fstream output(outputFile.c_str(), ios::out);
    if(!output)
    {
        cout<<"Can not open output file "<<outputFile<<endl;
      //  cerr << "Error: " << strerror(errno);
        exit(0);
    }
  //  srand(time(NULL));
    float r = 1.0;//0.7;
    float g = 1.0;//0.3;// 1;// 0.3;
    float b = 0.0;//1.0;
    
    
    if (outputFile.find("backtracked") != std::string::npos) {
        
        r = 1.0;//0.7;
        g = 0.2;//0.3;// 1;// 0.3;
        b = 0.2;//1.0;
        
        
        //Temp
        r = 0.0;//0.7;
        g = 1.0;//0.3;// 1;// 0.3;
        b = 0.7;//1.0;
        
        
    }
    else if (outputFile.find("smoothed") != std::string::npos) {
        r = 0.3;//0.7;
        g = 1.0;//1.0;//0.3;// 1;// 0.3;
        b = 0.0;//1.0;
        
        //Temp
        //   r = 0.5;//0.7;
        //  g = 0.3;//1.0;//0.3;// 1;// 0.3;
        // b = 0.1;//1.0;
    }
    else if (outputFile.find("traced") != std::string::npos) {
        r = 0.6;//0.7;
        g = 0.5;//0.3;// 1;// 0.3;
        b = 1.0;//1.0;
    }
    else {
        
        r = 1.0;//0.7;
        g = 0.0;//0.3;// 1;// 0.3;
        b = 0.0;//1.0;
    }
    
    string  filamentHeader = filamentName + "_" + to_string(filamentId);
    
    output <<"<marker_set name=\"Filament_"<<filamentHeader<<"\">\n";
    
    float radius = 30.0;
    if (coordinates.size() < 2) {
        radius = 38.0;
    }
    
    for (int i = 0; i < (int)coordinates.size(); i++) {
        
        Coordinate coordinate = coordinates.at(i);
        
        output <<"<marker id=\""<<i<<"\""<<" x=\""<<coordinate.xCor<<"\" "<<"y=\""<<coordinate.yCor<<"\" "<<"z=\""<<coordinate.zCor<<"\" "<<"r=\""<<r<<"\" g=\""<<g<<"\" b=\""<<b<<"\" radius=\""<<radius<<"\"/>";
        output<<"\n";
    }
    
    for (int i = 0; i < (int)coordinates.size() - 1; i++) {
        output <<"<link id1=\""<<i<<"\""<<" id2=\""<<i + 1<<"\""<<" r=\""<<r<<"\" g=\""<<g<<"\" b=\""<<b<<"\" radius=\"30.0\"/>";
        output<<"\n";
    }
    
    output <<"</marker_set>";
    // printf("Writting data in CMM format complete\n");
}

void MRC::writeInPDBFormat(string outputFile, vector<Coordinate> coordinates)
{
    //ofstream file(fileName.c_str());
    fstream output(outputFile.c_str(), ios::out);
    if(!output)
    {
        cout<<"Can not open output file "<<outputFile<<endl;
       // cerr << "Error: " << strerror(errno);
       // exit(1);
    }
    
    for (int i = 0; i < coordinates.size(); i++) {
        
        Coordinate coordinate = coordinates.at(i);
        
        // int xCord = (coordinate.xCor * 100);
        // int yCord = (coordinate.yCor * 100);
        // int zCord = (coordinate.zCor * 100);
        
        
        // coordinate.xCor = xCord/100;
        // coordinate.yCor = yCord/100;
        // coordinate.zCor = zCord/100;
        
        
        
        output<< std::setw(4);
        output<< std::left;
        output <<"ATOM";
        
        output<<"  ";
        
        output<< std::setw(5);
        output<< std::right;
        output <<"235";
        
        
        output<< std::setw(4);
        output<< std::right;
        output <<"CB";
        
        output<<"  "; //1 more space due to left alighnment;
        
        
        //int dd = (int)(denseVoxel.density * 100);
        output<< std::setw(3);
        output<< std::right;
        output <<"THR";
        // output <<dd;
        
        output<<" ";
        
        output<<"A";
        
        //   output<< std::fixed << std::setprecision(2); // ??? Added later
        output<< std::setw(4);
        output<< std::right;
        output <<100;
        // output <<dd;
        
        output<<"    "; // 4 spaces before x starts
        
        //295  CB  THR A  35";
        //std::ofstream << std::fixed;
        output<< std::fixed << std::setprecision(2);
        output<< std::setw(8);
        output<< std::right;
        // cout<<denseVoxel.xValue<<"\n";
        output <<coordinate.xCor;
        
        output<< std::fixed << std::setprecision(2);
        output<< std::setw(8);
        output<< std::right;
        output<<coordinate.yCor;
        
        output<< std::fixed << std::setprecision(2);
        output<< std::setw(8);
        output<< std::right;
        output<<coordinate.zCor;
        
        output<<"  1.00 15.54           C\n";
        //output <<"ATOM   295  CB  THR A  35"<<denseVoxel.xValue<<" "<<denseVoxel.yValue<<" "<<denseVoxel.zValue<<"#1.00 15.54           C\n";
    }
    
    printf("Writting to PDB format is complete\n");
}

void MRC::printDensityAtIndex(const Index& index)
{
  int x = index.xIndex;
  int y = index.yIndex;
  int z = index.zIndex;

  float density = cube[x][y][z];
  /// string statement = "The density at point [" + xIndex + "][" + yIndex + "][" + zIndex + "] = ";
  cout << "The density at point [" << x << "][" << y << "][" << z << "] = " << density << endl;
}

void printDensityAtIndex(float density, int x, int y, int z)
{
  cout << "The density at point [" << x << "][" << y << "][" << z << "] = " << density << endl;
}

void printAllDensitiesInCube(MRC& mrc)
{
  cout << "****************************************" << endl;
  for(int i = 0; i < mrc.nz; i++)
    {
		for(int j = 0; j < mrc.ny; j++)
		{
			for(int k = 0; k < mrc.nx; k++)
			{
				float density = mrc.cube[k][j][i];
				printDensityAtIndex(density, k, j, i);
				// cin.get();
			}
		}
    }
  cout << "****************************************" << endl;
}

void printDensityAtIndexFromCoordinate(const Coordinate& coordinate, MRC& mrc)
{
	
}

float MRC::getCellA0() const
{
	return xLength;
}

float MRC::getCellA1() const
{
	return yLength;
}

float MRC::getCellA2() const
{
	return zLength;
}

void MRC::printVoxelSize() const
{
	float xVoxelSize = xLength / mx;
    float yVoxelSize = yLength / my;
    float zVoxelSize = zLength / mz;    
    cout << "Voxel Size: " << endl;
    cout << "CELLA[0]/mx = " << xLength << "/" << mx << " = " << xVoxelSize << " angstroms on X" << endl
         << "CELLA[1]/my = " << yLength << "/" << my << " = " << yVoxelSize << " angstroms on Y" << endl
         << "CELLA[2]/mz = " << zLength << "/" << mz << " = " << zVoxelSize << " angstroms on Z" << endl;
}

float MRC::getVoxelSize() const
{
	float xVoxelSize = xLength / mx;
	float yVoxelSize = yLength / my;
	float zVoxelSize = zLength / mz;
	if(xVoxelSize == yVoxelSize && xVoxelSize == zVoxelSize && yVoxelSize == zVoxelSize)
	{
		float generalVoxelSize = xVoxelSize;
		return generalVoxelSize;
	}
	
	return xVoxelSize;
}

int MRC::convertAngstromsToVoxels(float angstromDistance) const
{
	// Add a value of 0.5 for use in symmetric rounding when casting float to int.
  	float numVoxels = (angstromDistance / this->getVoxelSize());
	cout << "Floating-Point Real Value: " << numVoxels << endl;
	numVoxels += 0.5;
	numVoxels = static_cast<int>(numVoxels);
	return numVoxels;
}

#pragma mark - Destructor
MRC::~MRC()
{
  for(int i = 0; i < nx; i++)
    for(int j = 0; j < ny; j++)
      {
		delete [] cube[i][j];
      }

  for(int i = 0; i < nx; i++)
    delete [] cube[i];

  delete [] cube;
}

void MRC::readMRCandSeeds(string inputMrcFilePath, string inputCmmFilePath, vector<Coordinate>& seeds)
{
    this->readMRCFile(inputMrcFilePath);
    this->readCoordinateFromCMMFile(seeds, inputCmmFilePath);
    cout << "Size of Seeds: " << seeds.size() << endl;
    cout << "The number of columns along the fast x-axis of the 3D data array cube is " << this->nx << endl
         << "The number of rows along the medium-speed y-axis of the 3D data array cube is " << this->ny << endl
         << "The number of sections along the slow z-axis of the 3D data array cube is " << this->nz << endl;
    // printAllDensitiesInCube(mrc);
    
    // Practice with writing to files of type .pdb so as to gain experience necessary for CS595 Python Assignment.
    // writeInPDBFormat("seeds.pdb", seeds);
}

void MRC::convertCoordinatesToIndices(vector<Coordinate> seeds)
{
    // Practice with conversion of coordinates to cube indices and access of seed point density data of each filament.
    vector<Coordinate>::iterator sitr = seeds.begin();
    Index testIndex;
    while(sitr != seeds.end())
    {
        this->getIndexFromCoordinate((float)(*sitr).xCor, (*sitr).yCor, (float)(*sitr).zCor, testIndex);
        cout << "Coordinate Values: " << endl
             << "X: " << (*sitr).xCor << endl
             << "Y: " << (*sitr).yCor << endl
             << "Z: " << (*sitr).zCor << endl << endl;
        cout << "Corresponding Index Values: " << endl
             << "X: " << testIndex.xIndex << endl
             << "Y: " << testIndex.yIndex << endl
             << "Z: " << testIndex.zIndex << endl << endl;
        this->printDensityAtIndex(testIndex);
        
        sitr++;
        cin.get();
    }
}

int main(int argc, char** argv)
{
    string inputMrcFilePath = (std::string)argv[1];
    string inputCmmFilePath = (std::string)argv[2];
    vector<Coordinate> seeds;
    MRC mrc;
    
    // Read input MRC file and associated Coordinate seeds.
    mrc.readMRCandSeeds(inputMrcFilePath, inputCmmFilePath, seeds);

    // Practice with reading of, access of, printing of, and conversion of voxel values to angstrom units along each axis.
    mrc.printVoxelSize();
    float voxelSize = mrc.getVoxelSize();

    // Convert seed Coordinate values to corresponding Index values and print results.
    // mrc.convertCoordinatesToIndices(seeds);
    
    /**
    string filamentExample;
    cout << "Enter example filament name: ";
    cin >> filamentExample;
    cout << endl;
    
    mrc.readCoordinateFromCMMFile(filament, filamentExample);
    int c = 0;
    for(vector<Coordinate>::iterator fitr = filament.begin(); fitr != filament.end(); fitr++)
    {
        cout << (*fitr);
        c++;
    }
    */

    // Load all existing filaments into vector collection for pruning of unnecessary coordinates
    vector<Coordinate> filament;
    vector<vector<Coordinate>> filaments(seeds.size()); // One filament per starting seed point, so filaments gets size seeds.size()
    int i = 0;
    for(vector<vector<Coordinate>>::iterator fitr = filaments.begin(); fitr != filaments.end(); fitr++) 
    {
        string filamentName = "/home/jhessefo/git/MRC/Chimera_MRC_Filament_Tracer/Output/finals_filament_smoothed/filament_smoothed_" + to_string(i) + ".cmm";
        mrc.readCoordinateFromCMMFile(filaments.at(i), filamentName);
        
        cout << endl << "Grabbing coordinates of filament " << i << "..." << endl;
        for(vector<Coordinate>::iterator citr = filaments[i].begin(); citr != filaments[i].end(); citr++) // iterate through each filament in filaments and print coordinate of the filament
        	cout << (*citr);
    	i++;
    }

    // Test: Measure average density along a single filament within voxel radius and with each measure check to see if average density has dropped below mean
    double meanDensity = mrc.meanDensity();
    // meanDensity = abs(meanDensity);
    int numCoords = filaments.at(1).size();
    vector<Coordinate>::iterator citr = filaments.at(1).begin();
    for(int j = 0; j < numCoords - 1; j++)
      {
	cout << "Coordinate " << j << ": " << filaments.at(1).at(j)
	     << "Coordinate " << j+1 << ": " << filaments.at(1).at(j+1);
	float filamentDistance = calculateDistance(filaments.at(1).at(j), filaments.at(1).at(j+1));
	float densityRadius = filamentDistance / 2;
	int voxelRadius = mrc.convertAngstromsToVoxels(densityRadius);
	cout << "Distance between coordinate " << j << " and coordinate " << j+1  << ": " << filamentDistance << endl
	     << "Density Calculation Radius: " << densityRadius << " angstroms" << endl
             << "Voxel Radius Equivalent: " << voxelRadius << endl;

	Index testIndex;
	mrc.getIndexFromCoordinate((float)filaments.at(1).at(j).xCor, filaments.at(1).at(j).yCor,(float)filaments.at(1).at(j).zCor, testIndex); 
	cout << "Corresponding Cube Indices: " << testIndex;
	double totalDensity = 0;
        for(int z = 0; z < voxelRadius; z++)
	  {
	    for(int y = 0; y < voxelRadius; y++)
	      {
		for(int x = 0; x < voxelRadius; x++)
		  {
		    if(testIndex.xIndex + x < mrc.nx && testIndex.yIndex + y < mrc.ny && testIndex.zIndex + z < mrc.nz)
		      {
			totalDensity += mrc.cube[testIndex.xIndex + x][testIndex.yIndex + y][testIndex.zIndex + z];
			// cout << "totalDensity = " << totalDensity << endl;
		      }
		    //totalDensity += mrc.cube[testIndex.xIndex - x][testIndex.yIndex - y][testIndex.zIndex - z];
		 }
	      }
	  }
	float averageDensity = (float) totalDensity / (voxelRadius * voxelRadius * voxelRadius); // Dimensions multiplied then squared to account for multi-directional forward/backward accumulation of densities
	// averageDensity = abs(averageDensity);
	cout << "Average Density: " << averageDensity << endl << endl;
  if(averageDensity > meanDensity)
    {
      cout << "Well, it looks like " << averageDensity << " is greater than " << meanDensity << ", so it seems we can determine precise termination is here at Coordinate " << j << endl;
      filaments[1].erase(citr, filaments[1].end());
      cout << "Final marking point of this filament: Marker (" << j << "): " << (*citr);
      break;
    }
  citr++;
 }
	/**	
   	int numCoords = seeds.size();
	for(int i = 0; i < numCoords - 1; i++)
	{
		float filamentDistance = calculateDistance(seeds.at(i), seeds.at(i+1));
		float densityRadius = filamentDistance / 2;
		int voxelRadius = mrc.convertAngstromsToVoxels(densityRadius);
		cout << "Distance between coordinate " << i << " and coordinate " << i+1 << ": " << filamentDistance << endl
			 << "Density Calculation Radius: " << densityRadius << " angstroms" << endl
			 << "Voxel Radius Equivalent: " << voxelRadius << endl << endl;
	}
	*/

    /**
     * Thinking this through: What do I need to do from here?
     * First, what can I do already?
     * - Access density values at any point and index on cube
     * - Convert coordinate values to cube index equivalents
     * - Read filament data from marker cmm file
     * - Calculate the distance between two coordinates, derive desired radius value from this, and convert the radius distance in angstroms to the voxel equivalent integer
     * Second, what do I do from here?
     * - Use voxel radius equivalent as number of indices to move in all dimensions, x, y, and z, in order to create a "cube minor" from which to take the average of all densities in the radius
     * - Take average density of entire cube and use as minimal boundary value for the average densities of cube minors
     * - If average cube minor density < boundary -> set marker at current coordinate value to determine terminal point of current filament, and begin examination of following filament
     * - Read filaments from output folder by incrementing integer to append to file name string before opening filament file of current loop iteration
     * - Vector of vectors of coordinates to represent all traced filaments? vector<vector<Coordinate>> 
     */

    return 0;
}

/**
#pragma mark - ReadNTrace Filament
void MRC::readSeedPointsNTraceFilaments(string seedPointsFilename)
{
    vector<Coordinate> seedpoints;

    readCoordinateFromCMMFile(seedpoints, seedPointsFilename);
    sort(seedpoints.begin(), seedpoints.end(), sortBasedOnZCoord);

    numberOfFilament =   (int)seedpoints.size();

    int firstOnZ = 0;
    int lastOnZ;

    cout<<"Sorting Seed Points...."<<endl;
    cout<<"-------------------------------------------"<<endl;
    isFirstFilemantOnPlane[0] = true;
    for (int i = 1; i < numberOfFilament; i++) {
        isFirstFilemantOnPlane[i] = false;
        if (fabs(seedpoints.at(i).zCor -  seedpoints.at(i - 1).zCor) > getApixZ() * 8){
            lastOnZ = i - 1;
            cout<<"Seed Plane: "<<firstOnZ<< " " <<lastOnZ << endl;
            sort(seedpoints.begin() + firstOnZ, seedpoints.begin() + lastOnZ + 1, sortBasedOnXCoord);
            firstOnZ = i;
            isFirstFilemantOnPlane[i] = true;
        }
    }

    sort(seedpoints.begin() + firstOnZ, seedpoints.begin() + numberOfFilament, sortBasedOnXCoord);

     if (filamentDirection == Ascending) {
         mrcEndY =  ny - 1;//2000;
     }
     else {
         mrcEndY = 1;//  90;
     }
    cout<<"-------------------------------------------"<<endl;

    cout<<"\nTracing start........."<<endl;
    cout<<"Phase 1: tracing using only density"<<endl;

    traceFilamentsWithStepSize(seedpoints, 2);
    int overlappingSegments2 = getNumberOfOverlappingSegments();
    //printf("Overlapping Segments for shift of 2: %d\n", overlappingSegments2);
  // readAndCalculateAccuracy();

    traceFilamentsWithStepSize(seedpoints, 3);
    int overlappingSegments3 = getNumberOfOverlappingSegments();
   // printf("Overlapping Segments forshift of 3: %d\n", overlappingSegments3);
  //  readAndCalculateAccuracy();

    if (overlappingSegments2 < overlappingSegments3) {
       // printf("Selecting Best x shift");
        traceFilamentsWithStepSize(seedpoints, 2);
    }


   // readAndCalculateAccuracy();
    cout<<"Phase-1 completed"<<endl<<endl;
}
**/
