/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 3.1                               *
 *                                                                       *
 * "generateMassMatrix" utility , Copyright (C) 2007 CMU, 2009 MIT,      *
 *                                              2016 USC                 *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic                                           *
 * http://www.jernejbarbic.com/code                                      *
 *                                                                       *
 * Research: Jernej Barbic, Fun Shing Sin, Daniel Schroeder,             *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
 *          Zumberge Research and Innovation Fund at USC                 *
 *                                                                       *
 * This utility is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this utility in the file LICENSE.txt                    *
 *                                                                       *
 * This utility is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

/*
  Generates the mass matrix for a given volumetric mesh (cubic or tet mesh).
*/

#include <fstream>
#include "generateMassMatrix.h"
#include "StVKStiffnessMatrix.h"
#include "StVKElementABCDLoader.h"
#include "ARPACKSolver.h"
#include "tetMesher.h"

using namespace std;

void printHelp(){
	// TODO
	std::cout << "FEM Faust Physical Model Generator\n";
}

/*
TODO:
- list of args for doc (finish printHelp):
	1: .obj file: the model
- compute reflectance option
- set object name
- set min and max freqs
*/

int main(int argc, char ** argv)
{
	const char *modelFileName = ""; // .obj file name
  std::string objectName = "modalModel"; // generated object name
  double materialProperties[3] = {70E9, 0.35, 2700}; // young's modulus, poisson's ratio and density - default: aluminum
	std::vector<int> exPos;
  bool debugMode = false; // debug mode activated
	bool showFreqs = false; // hide or show frequencies in the selected range
	float modesMinFreq = 20; // lowest mode freq
  float modesMaxFreq = 10000; // highest mode freq
  int targetNModes = 20; // number of synthesized modes
	int femNModes = 100; // number of synthesized modes
	int nExPos = -1; // number of excitation positions (default is max)
	int modesSelMode = 0; // mode to select modes (linear/max gains/cricital bands)

	/////////////////////////////////////
	// PARSING ARGUMENTS
	/////////////////////////////////////
  int currentArg = 0;
  if(argc > 1){
  	while(currentArg <= (argc-1)){
  		if(strcmp(argv[currentArg],"--infile") == 0){
  			currentArg++;
  			modelFileName = argv[currentArg];
				if(strcmp(modelFileName,"") == 0){
		  		cout << "No .obj file provided!\n";
		  		return 0;
		  	}
  		}
  		else if(strcmp(argv[currentArg],"--material") == 0){
  			currentArg++;
  			for(int i=0; i<3; i++){
  				if(currentArg > (argc-1)){
  					cout << "Not enough parameter properties!\n";
  					return 0;
  				}
  				if(strtod(argv[currentArg],NULL) != 0){
  					materialProperties[i] = strtod(argv[currentArg],NULL);
  				}
  				else{
  					cout << "Wrong material parameter: " << argv[currentArg] << "\n";
  					return 0;
  				}
  				if(i<2) currentArg++;
  			}
  		}
			else if(strcmp(argv[currentArg],"--expos") == 0){
  			currentArg++;
				while(strtod(argv[currentArg],NULL) != 0){
					exPos.push_back(strtod(argv[currentArg],NULL));
					currentArg++;
				}
				currentArg--;
  		}
  		else if(strcmp(argv[currentArg],"--debug") == 0){
  			debugMode = true;
  		}
			else if(strcmp(argv[currentArg],"--showfreqs") == 0){
  			showFreqs = true;
  		}
  		else if(strcmp(argv[currentArg],"--name") == 0){
  			currentArg++;
  			if(currentArg > (argc-1)){
  				cout << "--name: expecting an argument\n";
  				return 0;
  			}
  			objectName = argv[currentArg];
  		}
  		else if(strcmp(argv[currentArg],"--minmode") == 0){
  			currentArg++;
  			if(currentArg > (argc-1)){
  				cout << "--minmode: expecting an argument\n";
  				return 0;
  			}
  			if(strtod(argv[currentArg],NULL) != 0){
  				modesMinFreq = strtof(argv[currentArg],NULL);
  			}
  			else{
  				cout << "Min mode is not a float\n";
  				return 0;
  			}
  		}
  		else if(strcmp(argv[currentArg],"--maxmode") == 0){
  			currentArg++;
  			if(currentArg > (argc-1)){
  				cout << "--maxmode: expecting an argument\n";
  				return 0;
  			}
  			if(strtod(argv[currentArg],NULL) != 0){
  				modesMaxFreq = strtof(argv[currentArg],NULL);
  			}
  			else{
  				cout << "Max mode is not a float\n";
  				return 0;
  			}
  		}
			else if(strcmp(argv[currentArg],"--lmexpos") == 0){
  			currentArg++;
  			if(currentArg > (argc-1)){
  				cout << "--lmexpos: expecting an argument\n";
  				return 0;
  			}
  			if(strtod(argv[currentArg],NULL) != 0){
  				nExPos = strtod(argv[currentArg],NULL);
  			}
  			else{
  				cout << "Excitation position limit is not an int\n";
  				return 0;
  			}
  		}
  		else if(strcmp(argv[currentArg],"--nsynthmodes") == 0){
  			currentArg++;
  			if(currentArg > (argc-1)){
  				cout << "--nsynthmodes: expecting an argument\n";
  				return 0;
  			}
  			if(strtod(argv[currentArg],NULL) != 0){
  				targetNModes = strtod(argv[currentArg],NULL);
  			}
  			else{
  				cout << "Number of synthesized modes is not an int\n";
  				return 0;
  			}
  		}
			else if(strcmp(argv[currentArg],"--nfemmodes") == 0){
  			currentArg++;
  			if(currentArg > (argc-1)){
  				cout << "--nfemmodes: expecting an argument\n";
  				return 0;
  			}
  			if(strtod(argv[currentArg],NULL) != 0){
  				femNModes = strtod(argv[currentArg],NULL);
  			}
  			else{
  				cout << "Number of FEM modes is not an int\n";
  				return 0;
  			}
  		}
			else if(strcmp(argv[currentArg],"--modesel") == 0){
  			currentArg++;
  			if(currentArg > (argc-1)){
  				cout << "--modesel: expecting an argument\n";
  				return 0;
  			}
  			if(strtod(argv[currentArg],NULL) != 0){
  				modesSelMode = strtod(argv[currentArg],NULL);
  			}
  			else{
  				cout << "Modes selection mode is not an int\n";
  				return 0;
  			}
  		}
  		currentArg++;
  	}
  }
  else{
  	cout << "Missing arguments!\n\n";
  	printHelp();
  	return 0;
  }

	/////////////////////////////////////
	// RETRIEVING MODEL
	/////////////////////////////////////

	// loading mesh file
  if(debugMode){
		cout << "Loading the mesh file\n";
	}
  ObjMesh *objMesh = new ObjMesh(modelFileName);

	// generating 3D volumetric mesh from 2D mesh
  if(debugMode){
  	cout << "\nGenerating a 3D mesh with the following properties\n";
  	cout << "Young's modulus: " << materialProperties[0] << "\n";
  	cout << "Poisson's ratio: " << materialProperties[1] << "\n";
  	cout << "Density: " << materialProperties[2] << "\n";
  }
	// TODO: no way to prevent printing here (yet)
  TetMesh *volumetricMesh = TetMesher().compute(objMesh);
  // setting mesh material properties
  volumetricMesh->setSingleMaterial(materialProperties[0],materialProperties[1],materialProperties[2]);

	// computing mas matrix
  if(debugMode){
		cout << "Creating and computing mass matrix\n";
	}
  SparseMatrix *massMatrix;
  GenerateMassMatrix::computeMassMatrix(volumetricMesh, &massMatrix, true);

  // computing stiffness matrix
  if(debugMode){
		cout << "Creating and computing stiffness matrix\n";
	}
	StVKElementABCD *precomputedIntegrals = StVKElementABCDLoader::load(volumetricMesh);
	StVKInternalForces *internalForces = new StVKInternalForces(volumetricMesh, precomputedIntegrals);
	SparseMatrix *stiffnessMatrix;
	StVKStiffnessMatrix *stiffnessMatrixClass = new StVKStiffnessMatrix(internalForces);
	stiffnessMatrixClass->GetStiffnessMatrixTopology(&stiffnessMatrix);
	double *zero = (double*) calloc(3 * volumetricMesh->getNumVertices(), sizeof(double));
	stiffnessMatrixClass->ComputeStiffnessMatrix(zero, stiffnessMatrix);

	/////////////////////////////////////
	// EIGEN ANALYSIS
	/////////////////////////////////////

	// temporary variables for analysis
	//int numModes = stiffnessMatrix->Getn()-1; // number of computed modes: always max
	//int numModes = 110;
	double *eigenValues = (double*) malloc (sizeof(double)*femNModes);
  double *eigenVectors = (double*) malloc (sizeof(double)*stiffnessMatrix->Getn()*femNModes);

	// solver parameters
  double sigma = -1.0;
  int numLinearSolverThreads = 1; // by default only one thread but could try more...
	// starting solver
	if(debugMode){
		cout << "\nStarting the eigen solver\n";
		cout << femNModes << " modes will be computed for the FEM analysis\n\n";
	}
	ARPACKSolver generalizedEigenvalueProblem;
  int nconv = generalizedEigenvalueProblem.SolveGenEigShInv(stiffnessMatrix, massMatrix, femNModes, eigenValues, eigenVectors, sigma, numLinearSolverThreads,0);
  //int nconv = generalizedEigenvalueProblem.SolveGenEigReg(stiffnessMatrix, massMatrix, femNModes, eigenValues, eigenVectors, "LM", numLinearSolverThreads,0);

  if(nconv == femNModes){ // if analysis was successful...

		/////////////////////////////////////
		// COMPUTING MODE FREQS
		/////////////////////////////////////

		if(debugMode){
			printf("Computing modes frequencies\n\n");
		}
  	float modesFreqs[femNModes]; // modes freqs
		int lowestModeIndex = 0;
		int highestModeIndex = 0;
  	for(int i=0; i<femNModes; i++){
  		if (eigenValues[i] <= 0){
		  	modesFreqs[femNModes] = 0.0;
	  	}
	  	else{
		  	modesFreqs[i] = sqrt((float)eigenValues[i]) / (2 * M_PI);
	  	}
			if(modesFreqs[i] < modesMinFreq){
				lowestModeIndex++;
			}
			if(modesFreqs[i] < modesMaxFreq && (highestModeIndex-lowestModeIndex) <= targetNModes && highestModeIndex < volumetricMesh->getNumVertices()){
				highestModeIndex++;
			}
  	}

		// adjusting number of target modes to modes range
		int modesRange = highestModeIndex-lowestModeIndex;
		if(modesRange < targetNModes){
			targetNModes = modesRange;
		}

		// diplaying mode frequencies
		if(showFreqs){
			cout << "Mode frequencies between " << modesMinFreq << "Hz and " << modesMaxFreq << "Hz:\n";
			for(int i=0; i<targetNModes; i++){
				cout << modesFreqs[i+lowestModeIndex] << "\n";
			}
			cout << "\n";
		}

		/////////////////////////////////////
		// COMPUTING GAINS
		/////////////////////////////////////

		if(debugMode){
			cout << "Computing modes gains for modes between " << modesMinFreq << "Hz and " << modesMaxFreq << "Hz\n\n";
		}
		if(exPos.size() > 0){ // if exPos specified, then retrieve number of ex positions
			nExPos = exPos.size();
		}
		if(nExPos == -1 || nExPos > volumetricMesh->getNumVertices()){ // if nExPos not specified, then max number of exPos
			nExPos = volumetricMesh->getNumVertices();
		}
		int exPosStep = volumetricMesh->getNumVertices()/nExPos; // to skip excitation positions
		float modesGains[nExPos][targetNModes]; // modes gains matrix
		for(int i=0; i<nExPos; i++){ // i = excitation position
			float maxGain = 0; // for normalization
			for(int j=0; j<targetNModes; j++){ // j = current mode
				if(exPos.size()>0){ // if expos was defined, then retrieve data
					modesGains[i][j] = sqrt(
						pow(eigenVectors[(j+lowestModeIndex)*stiffnessMatrix->Getn()+(exPos[i]*3)],2) +
						pow(eigenVectors[(j+lowestModeIndex)*stiffnessMatrix->Getn()+(exPos[i]*3)+1],2) +
						pow(eigenVectors[(j+lowestModeIndex)*stiffnessMatrix->Getn()+(exPos[i]*3)+2],2)
					);
				}
				else{ // otherwise choose linear ex pos
					modesGains[i][j] = sqrt(
						pow(eigenVectors[(j+lowestModeIndex)*stiffnessMatrix->Getn()+(i*exPosStep*3)],2) +
						pow(eigenVectors[(j+lowestModeIndex)*stiffnessMatrix->Getn()+(i*exPosStep*3)+1],2) +
						pow(eigenVectors[(j+lowestModeIndex)*stiffnessMatrix->Getn()+(i*exPosStep*3)+2],2)
					);
				}
				if(modesGains[i][j]>maxGain) maxGain = modesGains[i][j]; // saving max gain for normalization
			}
			// normalizing gains for current position
			for(int j=0; j<targetNModes; j++){
				modesGains[i][j] = modesGains[i][j]/maxGain;
				/*
				if(i==0){
					cout << modesFreqs[lowestModeIndex+j] << ":\t" << modesGains[i][j] << "\n";
				}
				*/
			}
		}

		// TODO missing section here for mode selection. In both cases, matrix of indices
		// will have to be computed and can use the same variable. This will have to be
		// taken into account when generating the Faut files.
		// TODO: also need to display the result of this operation

		/////////////////////////////////////
		// GENERATING FAUST FILE
		/////////////////////////////////////

		// TODO: say something about the model that will be generated (parameters available, etc.)

		std::string faustFileName;
		faustFileName.append(objectName).append(".lib");
		std::ofstream faustFile(faustFileName.c_str());
		faustFile << "import(\"stdfaust.lib\");\nimport(\"lpm.lib\");\n\n";
		faustFile << "nModes = " << targetNModes << ";\n";
		if(nExPos > 1) faustFile << "nExPos = " << nExPos << ";\n";

		// TODO generating static model for now
		faustFile << "modesFreqs(n) = ba.take(n+1,(";
		for(int i=0; i<targetNModes; i++){
			faustFile << modesFreqs[lowestModeIndex+i];
			if(i<(targetNModes-1)) faustFile << ",";
		}
		faustFile << "));\n";
		faustFile << "modesGains(p,n) = waveform{";
		for(int i=0; i<nExPos; i++){
			for(int j=0; j<targetNModes; j++){
				faustFile << modesGains[i][j];
				if(j<(targetNModes-1)) faustFile << ",";
			}
			if(i<(nExPos-1)) faustFile << ",";
		}
		faustFile << "},int(p*nModes+n) : rdtable;\n\n";
		faustFile << "modalModel(exPos,t60,t60DecayRatio,t60DecaySlope) = _ <: par(i,nModes,*(modesGains(int(exPos),i)) : modeFilter(modesFreqs(i),modesT60s(i))) :> /(nModes)\n";
		faustFile << "with{\n";
		float freqDiff = modesFreqs[lowestModeIndex]/modesFreqs[highestModeIndex];
		faustFile << "modesT60s(i) = t60*pow(1-(modesFreqs(i)/" << modesFreqs[highestModeIndex] << " - " << freqDiff << ")*(t60DecayRatio + " << freqDiff << "),t60DecaySlope);\n";
		faustFile << "};\n";
		faustFile.close();

		/*
		faustFile << "modesFreqs(n) = ba.take(n+1,(";
		if(exPos == -1){
			for(int j=0; j<targetNModes; j++){
				//faustFile << modesFreqs[finalModesI[i*exStep][j]+lowestModeIndex];
				faustFile << modesFreqs[lowestModeIndex+j];
				if(j<(targetNModes-1)) faustFile << ",";
			}
		}
		else{
			for(int j=0; j<targetNModes; j++){
				//faustFile << modesFreqs[finalModesI[exPos][j]+lowestModeIndex];
				if(j<(targetNModes-1)) faustFile << ",";
			}
		}
		faustFile << "));\n";
		faustFile << "modesGains(p,n) = waveform{";
		if(exPos != -1){
			for(int j=0; j<targetNModes; j++){
				//faustFile << modesGainsMatrix[finalModesI[exPos][j]][exPos];
				if(j<(targetNModes-1)) faustFile << ",";
			}
		}
		else{
			for(int i=0; i<nExPos; i++){
				for(int j=0; j<targetNModes; j++){
					//faustFile << modesGainsMatrix[finalModesI[i*exStep][j]][i*exStep]; // 209 & 629
					faustFile << modesGains[i*exStep][j];
					if(j<(targetNModes-1)) faustFile << ",";
				}
				if(i<(nExPos-1)) faustFile << ",";
			}
		}
		faustFile << "},int(p*nModes+n) : rdtable;\n\n";
		//faustFile << "modalModel(scaleFreq,exPos,exSpread,maxT60,shortestT60) = _ <: par(i,nModes,*(mgs(i)) : modeFilter(modesFreqs(pos,i),modesT60(i))) :> /(nModes)\n";
		faustFile << "modalModel(scaleFreq,exPos,exSpread,maxT60,t60Decay,t60Slope) = _ <: par(i,nModes,*(modesGains(pos,i)) : modeFilter(modesFreqs(i),modesT60(i))) :> /(nModes)\n";
		faustFile << "with{\n";
		faustFile << "pos = int(exPos*nPosEx);\n";
		//faustFile << "modesT60(i) = abs(modesGains(pos,i))*maxT60*(1-((i+1)/nModes)*shortestT60);\n";
		faustFile << "modesT60(i) = maxT60*pow(1-modesFreqs(i)/" << modesMaxFreq << "*t60Decay,t60Slope);\n";
		//faustFile << "mgs(i) = (1-(modesFreqs(pos,i)/" << modesMaxFreq << ")*shortestT60);\n";
		faustFile << "};\n";
		*/
		//faustFile.close();

			/*
			int ex = 1;
			for(int i=0; i<numModes; i++){
				std::cout << modesFreqs[i] << "\t" << sqrt(pow(modesTemp[i*stiffnessMatrix->Getn()+(ex*3)],2) + pow(modesTemp[i*stiffnessMatrix->Getn()+(ex*3)+1],2) + pow(modesTemp[i*stiffnessMatrix->Getn()+(ex*3)+2],2)) << "\n";
				//std::cout << modesFreqs[i] << "\t" << sqrt(pow(modesTemp[ex*numModes+(i*3)],2) + pow(modesTemp[ex*numModes+(i*3)+1],2) + pow(modesTemp[ex*numModes+(i*3)+2],2)) << "\n";
			}
			*/

			/*
			for(int i=0; i<numModes; i++){
				for(int j=0; j<stiffnessMatrix->Getn(); j++){
					//std::cout << modesTemp[i*stiffnessMatrix->Getn()+j] << ",";
					std::cout << modesTemp[j*numModes+i] << ",";
				}
				std::cout << "\n";
			}
			*/

			// vertices / rLin / modes
			//ModalMatrix *modalMatrix = new ModalMatrix(stiffnessMatrix->Getn()-100, numModes-100, modesTemp);


/*
  		if(debugMode) printf("\nComputing modes amplitude matrix\n");
  		Eigen::MatrixXd modesMatrix(numModes,numModes);
  		for(int i=0; i<numModes; i++){
	  		for(int j=0; j<numModes; j++){
		  		modesMatrix(i,j) = modesTemp[(stiffnessMatrix->Getn()*j)+i];
	  		}
  		}
  		Eigen::MatrixXd modesMatrixInverse(numModes,numModes);
  		//modesMatrixInverse = modesMatrix.inverse();
			modesMatrixInverse = modesMatrix;

  		// detecting lowest and highest mode index
  		int lowestModeIndex = 0;
  		while(modesFreqs[lowestModeIndex] <= modesMinFreq){
  			lowestModeIndex++;
  		}
  		int highestModeIndex = 0;
  		while(modesFreqs[highestModeIndex] <= modesMaxFreq){
  			highestModeIndex++;
  		}
			int nModesRange = highestModeIndex-lowestModeIndex;

			// in case there are less modes in the defined range than desired
			if(targetNModes>nModesRange){
				targetNModes = nModesRange;
			}

  		if(debugMode){
  			std::cout << "\nGenerating a faust object with the following constraints:\n";
  			std::cout << "Lowest mode frequency: " << modesMinFreq << "\n";
  			std::cout << "Highest mode frequency: " << modesMaxFreq << "\n";
  			std::cout << "Number of available modes in this range: " << nModesRange << "\n";
  			std::cout << "Number of target modes: " << targetNModes << "\n";
  		}

  		// calculating critical bands
  		float modesRange = modesMaxFreq - modesMinFreq;
  		float criticalBands[targetNModes+1];
  		for(int i=0; i<=targetNModes; i++){
  			criticalBands[i] = pow(i/(float)targetNModes,2)*modesRange+modesMinFreq;
  		}

  		// ordering modes indices in vectors by critical bands
  		std::vector<int> modesIBands[targetNModes];
  		int cnt = 0;
  		int maxRowSize = 0;
  		for(int i=0; i<targetNModes; i++){
  			while(modesFreqs[cnt] < criticalBands[i+1] && cnt < numModes){
  				if(modesFreqs[cnt] >= criticalBands[i]){
  					modesIBands[i].push_back(cnt-lowestModeIndex);
  				}
  				cnt++;
  			}
  			if(modesIBands[i].size()>maxRowSize){
  				maxRowSize = modesIBands[i].size();
  			}
  		}

  		// previous vector array is turned into a matrix for export
  		int iMatrix[targetNModes][maxRowSize];
  		for(int i=0; i<targetNModes; i++){
  			for(int j=0; j<maxRowSize; j++){
  				if(modesIBands[i].size() == 0 || j >= modesIBands[i].size()){
  					iMatrix[i][j] = -1;
  				}
  				else{
  					iMatrix[i][j] = modesIBands[i][j];
  				}
  			}
  		}

			int n2dVert = stiffnessMatrix->Getn()/3;
			float modesGainsMatrix[nModesRange][n2dVert];
			float maxGains[n2dVert];
			for(int i=0; i<n2dVert; i++){
				for(int j=0; j<nModesRange; j++){
					//modesGainsMatrix[j][i] = modesMatrixInverse(lowestModeIndex+j,i) + modesMatrixInverse(lowestModeIndex+j,i+1) + modesMatrixInverse(lowestModeIndex+j,i+2);
					//modesGainsMatrix[j][i] = sqrt(pow(modesMatrix(lowestModeIndex+j,i),2) + pow(modesMatrix(lowestModeIndex+j,i+1),2) + pow(modesMatrix(lowestModeIndex+j,i+2),2));
					//modesGainsMatrix[j][i] = sqrt(pow(modesMatrix(i,lowestModeIndex+j),2) + pow(modesMatrix(i+1,lowestModeIndex+j),2) + pow(modesMatrix(i+2,lowestModeIndex+j),2));
					//modesGainsMatrix[j][i] = sqrt(pow(modesTemp[(j+lowestModeIndex)*stiffnessMatrix->Getn()+(i*3)],2) + pow(modesTemp[(j+lowestModeIndex)*stiffnessMatrix->Getn()+(i*3)+1],2) + pow(modesTemp[(j+lowestModeIndex)*stiffnessMatrix->Getn()+(i*3)+2],2));
					modesGainsMatrix[j][i] = sqrt(pow(modesTemp[i*numModes+((j+lowestModeIndex)*3)],2) + pow(modesTemp[i*numModes+((j+lowestModeIndex)*3)+1],2) + pow(modesTemp[i*numModes+((j+lowestModeIndex)*3)+2],2));
					float currentABSGain = std::abs(modesGainsMatrix[j][i]);
					if(maxGains[i] < currentABSGain){
						maxGains[i] = currentABSGain;
					}
				}
			}
			// normalizing modes gains
			for(int i=0; i<n2dVert; i++){
				for(int j=0; j<nModesRange; j++){
					modesGainsMatrix[j][i] = modesGainsMatrix[j][i]/maxGains[i];
				}
			}

			// building the gain matrix
			int finalModesI[n2dVert][targetNModes];
			for(int h=0; h<n2dVert; h++){
  			std::vector<float> modesGainsBands[targetNModes];
  			for(int i=0; i<targetNModes; i++){
  				for(int j=0; j<modesIBands[i].size(); j++){
						modesGainsBands[i].push_back(fabs(modesGainsMatrix[modesIBands[i][j]][h]));
  				}
  				std::sort(modesGainsBands[i].begin(),modesGainsBands[i].end(), std::greater<float>());
  			}
	  		cnt = 0;
	  		for(int j=0; (j<maxRowSize && cnt<targetNModes); j++){
	  			for(int i=0; (i<targetNModes && cnt<targetNModes); i++){
	  				if(modesGainsBands[i].size() > j){
	  					int k = 0;
	  					while(modesGainsBands[i][j] != (float)modesGainsMatrix[modesIBands[i][k]][h] &&
	  						modesGainsBands[i][j] != (float)-modesGainsMatrix[modesIBands[i][k]][h] &&
	  						k < modesGainsBands[i].size()){
	  						k++;
	  					}
	  					finalModesI[h][cnt] = modesIBands[i][k];
	  					cnt++;
	  				}
	  			}
	  		}
			}

			if(exPosLimit == -1 || exPosLimit > n2dVert){
				exPosLimit = n2dVert;
			}
			int exStep = n2dVert/exPosLimit;
			if(debugMode){
				std::cout << "Number of excitation positions: " << exPosLimit << "\n" ;
			}

			// Select by amplitude
			/*
			std::vector<float> mg[n2dVert];
			std::vector<float> mgs[n2dVert];
			int mi[n2dVert][targetNModes];
			for(int i=0; i<n2dVert; i++){
				float maxG = 0;
			  for(int j=0; j<std::min(nModesRange,volumetricMesh->getNumVertices()); j++){
					std::cout << j << " ";
			    mg[i].push_back(sqrt(pow(modesTemp[(j+lowestModeIndex)*stiffnessMatrix->Getn()+(i*3)],2) + pow(modesTemp[(j+lowestModeIndex)*stiffnessMatrix->Getn()+(i*3)+1],2) + pow(modesTemp[(j+lowestModeIndex)*stiffnessMatrix->Getn()+(i*3)+2],2)));
					if(mg[i][j]>maxG) maxG = mg[i][j];
			  }
				for(int j=0; j<std::min(nModesRange,volumetricMesh->getNumVertices()); j++){
					mg[i][j] = mg[i][j]/maxG;
					mgs[i].push_back(mg[i][j]);
				}
			  std::sort(mgs[i].begin(),mgs[i].end(), std::greater<float>());
			  for(int j=0; j<targetNModes; j++){
			    mi[i][j] = 0;
			    while(mgs[i][j] != mg[i][mi[i][j]]){
			      mi[i][j]++;
			    }
			  }
			}


			for(int j=0; j<numModes; j++){
			  std::cout << j << ": " << mg[0][j] << "\n";
			}

			for(int j=0; j<targetNModes; j++){
			  std::cout << modesFreqs[mi[0][j]] << "\n";
			}
*/

// HERE IN
/*
if(nExPos == -1 || nExPos > volumetricMesh->getNumVertices()){
	nExPos = volumetricMesh->getNumVertices();
}
int exStep = volumetricMesh->getNumVertices()/nExPos;

			std::string faustFileName;
  		faustFileName.append(objectName).append(".lib");
			std::ofstream faustFile(faustFileName.c_str());
			faustFile << "import(\"stdfaust.lib\");\nimport(\"lpm.lib\");\n\n";
			faustFile << "nModes = " << targetNModes << ";\n";
			if(exPos != -1) faustFile << "nPosEx = 0;\n";
			else faustFile << "nPosEx = " << nExPos << ";\n";
			faustFile << "modesFreqs(n) = ba.take(n+1,(";
			if(exPos == -1){
				for(int j=0; j<targetNModes; j++){
					//faustFile << modesFreqs[finalModesI[i*exStep][j]+lowestModeIndex];
					faustFile << modesFreqs[lowestModeIndex+j];
					if(j<(targetNModes-1)) faustFile << ",";
				}
			}
			else{
				for(int j=0; j<targetNModes; j++){
					//faustFile << modesFreqs[finalModesI[exPos][j]+lowestModeIndex];
					if(j<(targetNModes-1)) faustFile << ",";
				}
			}
			faustFile << "));\n";
			faustFile << "modesGains(p,n) = waveform{";
			if(exPos != -1){
				for(int j=0; j<targetNModes; j++){
					//faustFile << modesGainsMatrix[finalModesI[exPos][j]][exPos];
					if(j<(targetNModes-1)) faustFile << ",";
				}
			}
			else{
				for(int i=0; i<nExPos; i++){
					for(int j=0; j<targetNModes; j++){
						//faustFile << modesGainsMatrix[finalModesI[i*exStep][j]][i*exStep]; // 209 & 629
						faustFile << modesGains[i*exStep][j];
						if(j<(targetNModes-1)) faustFile << ",";
					}
	  			if(i<(nExPos-1)) faustFile << ",";
				}
			}
			faustFile << "},int(p*nModes+n) : rdtable;\n\n";
			//faustFile << "modalModel(scaleFreq,exPos,exSpread,maxT60,shortestT60) = _ <: par(i,nModes,*(mgs(i)) : modeFilter(modesFreqs(pos,i),modesT60(i))) :> /(nModes)\n";
			faustFile << "modalModel(scaleFreq,exPos,exSpread,maxT60,t60Decay,t60Slope) = _ <: par(i,nModes,*(modesGains(pos,i)) : modeFilter(modesFreqs(i),modesT60(i))) :> /(nModes)\n";
			faustFile << "with{\n";
			faustFile << "pos = int(exPos*nPosEx);\n";
			//faustFile << "modesT60(i) = abs(modesGains(pos,i))*maxT60*(1-((i+1)/nModes)*shortestT60);\n";
			faustFile << "modesT60(i) = maxT60*pow(1-modesFreqs(i)/" << modesMaxFreq << "*t60Decay,t60Slope);\n";
			//faustFile << "mgs(i) = (1-(modesFreqs(pos,i)/" << modesMaxFreq << ")*shortestT60);\n";
			faustFile << "};\n";
			faustFile.close();
*/
// HERE OUT
  	}

  	// cleaning
  	free(eigenValues);
  	free(eigenVectors);
  	free(zero);
  	delete stiffnessMatrixClass;
  	stiffnessMatrixClass = NULL;
  	delete stiffnessMatrix;
  	stiffnessMatrix = NULL;
  	delete internalForces;
	internalForces = NULL;
	delete precomputedIntegrals;
  	precomputedIntegrals = NULL;
  	delete massMatrix;
  	massMatrix = NULL;
  	delete volumetricMesh;
  	volumetricMesh = NULL;
  	delete objMesh;
  	objMesh = NULL;


  	return 0;

/*







  double **modes_;
  int n3 = 3 * objMesh->getNumVertices();
  int numConstrainedVertices = 0;
  int * constrainedDOFs = (int*) malloc (sizeof(int) * 3 * numConstrainedVertices);
  *modes_ = (double*) calloc (numDesiredModes * n3, sizeof(double));
  *?


  for(int i=0; i<numDesiredModes; i++){
	  // insert zero rows into the computed modes
	  int oneIndexed = 1;
	  InsertRows(n3, &modesTemp[numRetainedDOFs*i], &((*modes_)[n3*i]),
				 3 * numConstrainedVertices, constrainedDOFs, oneIndexed);
  }

  int n3 = 3 * volumetricMesh->getNumVertices();
  double **frequencies_;
  double **modes_;
  *frequencies_ = (double*) calloc ((numDesiredModes+1), sizeof(double));
  *modes_ = (double*) calloc ((numDesiredModes+1) * n3, sizeof(double));


  for(int i=0; i<numDesiredModes*numRetainedDOFs; i++){
	  printf("Mode %d: %f\n",i,modesTemp[i]);
  }

  //double modesGains[numDesiredModes];
  //double modesFreqs[numDesiredModes];


  float **modesMatrix;
  float **invModesMatrix;

  Eigen::MatrixXd m(numDesiredModes,numDesiredModes);
  Eigen::VectorXd v(numDesiredModes);

  modesMatrix = new float*[numDesiredModes];
  invModesMatrix = new float*[numDesiredModes];

  for(int i=0; i<numDesiredModes; i++){
	  modesMatrix[i] = new float[numDesiredModes];
	  invModesMatrix[i] = new float[numDesiredModes];
  }

  for(int i=0; i<numDesiredModes; i++){
	  if(i==300) v(i) = 1;
	  else v(i) = 0;
	  for(int j=0; j<numRetainedDOFs-1; j++){
		  m(j,i) = modesTemp[(numRetainedDOFs*i)+j];
		  //modesMatrix[j][i] = modesTemp[(numRetainedDOFs*i)+j];
		  //printf("Row %d \t Col %d \t %f\n",i,j,modesTemp[(numRetainedDOFs*i)+j]);
	  }
  }

  Eigen::MatrixXd mInv;
  mInv = m.inverse();
  std::cout << mInv*v;
  //m.diagonal();


  for(int i=0; i<numDesiredModes; i++){
	  modesGains[i] = modesTemp[numRetainedDOFs*i];
	  if (frequenciesTemp[i] <= 0){
		  modesFreqs[numDesiredModes] = 0.0;
	  }
	  else{
		  modesFreqs[i] = sqrt((frequenciesTemp)[i]) / (2 * M_PI);
	  }
	  printf("Modes %d \t Freq: %f \t Gain: %f\n",i,modesFreqs[i],modesGains[i]);
  }

  printf("Freqs:\n");
  for(int i=0; i<numDesiredModes; i++){
	  printf("%f\n",modesFreqs[i]);
  }

  printf("Gains:\n");
  for(int i=0; i<numDesiredModes; i++){
	  printf("%f\n",modesGains[i]);
  }

  free(modesTemp);
  free(frequenciesTemp);
  */
}
