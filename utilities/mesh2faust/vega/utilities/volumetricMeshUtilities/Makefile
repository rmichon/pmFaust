ifndef VOLUTILS
VOLUTILS=VOLUTILS

ifndef CLEANFOLDER
CLEANFOLDER=VOLUTILS
endif

include ../../Makefile-headers/Makefile-header
R ?= ../..


# the object files to be compiled for this utility
VOLUTILS_OBJECTS=

# the libraries this utility depends on
#VOLUTILS_LIBS=volumetricMesh reducedStvk stvk reducedElasticForceModel reducedForceModel forceModel renderVolumetricMesh objMesh imageIO modalMatrix matrix matrixIO getopts loadList sparseMatrix sparseSolver sparseMatrix graph minivector getopts insertRows loadList camera openGLHelper

#VOLUTILS_LIBS=mesher windingNumber volumetricMesh stvk reducedElasticForceModel reducedForceModel forceModel renderVolumetricMesh sparseSolver sparseMatrix objMesh imageIO modalMatrix matrix matrixIO getopts insertRows loadList camera minivector

VOLUTILS_LIBS=mesher windingNumber volumetricMesh stvk sparseSolver sparseMatrix objMesh matrix matrixIO

# the headers in this utility
VOLUTILS_HEADERS=

# TODO: this ight need some cleaning
VOLUTILS_LINK=$(addprefix -l, $(VOLUTILS_LIBS)) $(ARPACK_LIB) $(SPOOLES_LIB) $(BLASLAPACK_LIB) $(PARDISO_LIB) $(STANDARD_LIBS)

VOLUTILS_OBJECTS_FILENAMES=$(addprefix $(R)/utilities/volumetricMeshUtilities/, $(VOLUTILS_OBJECTS))
VOLUTILS_HEADER_FILENAMES=$(addprefix $(R)/utilities/volumetricMeshUtilities/, $(VOLUTILS_HEADERS))
VOLUTILS_LIB_MAKEFILES=$(call GET_LIB_MAKEFILES, $(VOLUTILS_LIBS))
VOLUTILS_LIB_FILENAMES=$(call GET_LIB_FILENAMES, $(VOLUTILS_LIBS))

include $(VOLUTILS_LIB_MAKEFILES)

all: $(R)/utilities/volumetricMeshUtilities/generateMassMatrix

$(R)/utilities/volumetricMeshUtilities/generateMassMatrix: $(R)/utilities/volumetricMeshUtilities/generateMassMatrix.cpp
	$(CXXLD) $(LDFLAGS) $(INCLUDE) $(VOLUTILS_OBJECTS) $(R)/utilities/volumetricMeshUtilities/generateMassMatrix.cpp $(VOLUTILS_LINK) -lgeometryQuery -lmeshKey -lrational -o $@;

$(VOLUTILS_OBJECTS_FILENAMES): %.o: %.cpp $(VOLUTILS_LIB_FILENAMES) $(VOLUTILS_HEADER_FILENAMES)
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) $< -o $@

ifeq ($(CLEANFOLDER), VOLUTILS)
clean: cleanvolumetricMeshUtilities
endif

deepclean: cleanvolumetricMeshUtilities

cleanvolumetricMeshUtilities:
	$(RM) $(VOLUTILS_OBJECTS_FILENAMES) $(R)/utilities/volumetricMeshUtilities/generateMassMatrix

endif
