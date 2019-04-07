/**
 *	Implementation file to support Chimera visualization tool plugin.
 *	Allows the user to work with functions to deal with membrane and filaments of hair cell stereocilia MRC images.
 */

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "MRC_Test.h"

using namespace std;

#pragma mark - Tracing
void MRC::traceFilament(Index startVoxel, Index endVoxel, vector<Coordinate> &maxDensityVoxels, int maximumXShiftBoundingBox, int filamentID)
{
	Index newRefVoxelSize20;
	float avgDensitySize20;

	Index refVoxel = startVoxel;
}

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
        stringstream sStream(line); // Initialize a stringstream with the contents of an entire text line as stored in var line.
        // std::smatch match;
        
        vector <std::string> splitedStrings; // Create a vector of streams to store the line string stream when separated into different segments.
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

/**
 * getIndexFromCoordinate() to convert a physical coordinate's x, y, z values into corresponding x, y, and z indices for the cube 3D array/triple pointer variable.
 * Once we have the converted indices x, y, and z, we can access the densities at a certain point through usage of cube[x][y][z].
 */

void MRC::getIndexFromCoordinate(float xCoord, int yCoord, float zCoord, Index &seedIndex)
{
    seedIndex.xIndex = ((xCoord - xOrigin)/ getApixX()) - nxStart ;
    seedIndex.yIndex = ((yCoord - yOrigin)/ getApixY()) - nyStart ;
    seedIndex.zIndex = ((zCoord - zOrigin)/ getApixZ()) - nzStart ;
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

int main(int argc, char* args[])
{
	inputMrcFilePath = (std::string)argv[1]; // \\undergrad-cifs\Undergrad\jhessefo\win_user_profile\Downloads\All Averaged + Orginial Map-20190131T140815Z-001\Averaged_Original\Dataset__1__ZM_Box25-4_110613_H10_16_Inverted_Cropped_Rotated_ActinsAligned(Voxel size=9.47).mrc
	inputCmmFilePath = (std::string)argv[2]; // \\undergrad-cifs\Undergrad\jhessefo\win_user_profile\Downloads\All Averaged + Orginial Map-20190131T140815Z-001\Averaged_Original\H10_16_ShaftDensitiesAtPlane1190.cmm
	vector<Coordinate> seeds;
	MRC mrc;	
	
	mrc.readMRCFile(inputMrcFilePath);
	mrc.readCoordinateFromCMMFile(seeds, inputCmmFilePath);
	cout << "Size of Seeds: " << seeds.size() << endl;


	return 0;
}

void printDensityAtIndex(const Index& index)
{
	int xIndex = index.x;
	int yIndex = index.y;
	int zIndex = index.z;

	float density = cube[xIndex][yIndex][zIndex];
	/// string statement = "The density at point [" + xIndex + "][" + yIndex + "][" + zIndex + "] = ";
	cout << "The density at point [" << xIndex << "][" << yIndex << "][" << zIndex << "] = " << density << endl;
}