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
			if(modesFreqs[i] < modesMaxFreq && (highestModeIndex-lowestModeIndex) < targetNModes && highestModeIndex < volumetricMesh->getNumVertices()){
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
		faustFile << objectName << "(exPos,t60,t60DecayRatio,t60DecaySlope) = _ <: par(i,nModes,*(modesGains(int(exPos),i)) : modeFilter(modesFreqs(i),modesT60s(i))) :> /(nModes)\n";
		faustFile << "with{\n";
		float freqDiff = modesFreqs[lowestModeIndex]/modesFreqs[highestModeIndex];
		faustFile << "modesT60s(i) = t60*pow(1-(modesFreqs(i)/" << modesFreqs[highestModeIndex] << " - " << freqDiff << ")*(t60DecayRatio + " << freqDiff << "),t60DecaySlope);\n";
		faustFile << "};\n";
		faustFile.close();
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
}
