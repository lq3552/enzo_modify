/***********************************************************************
/
/  INITIALIZE A GALAXY SIMULATION
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:  Elizabeth Tasker, March 2004
/  major change: Qi Li, May 2016
/  minor change1: Qi Li, Sep 2016
/  minor change2: Qi Li, Oct 2016
/  PURPOSE:
/
/    Set up an isolated disk galaxy
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
//

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"

void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);


int MyGalaxyInitialize(FILE *fptr, FILE *Outfptr, 
			  HierarchyEntry &TopGrid, TopGridData &MetaData)
{  /* need change when using introducing chemestry! */
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *G0Name = "G0";
  char *ZetaName = "Zeta";
  char *B1Name = "B1";
  char *B2Name = "B2";
  char *B3Name = "B3";
  char *PhiName = "Phi";
  char *MetalName = "Metal_Density";
  char *MetalIaName = "MetalSNIa_Density";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, disk, i,NumberOfSubgridZones[MAX_DIMENSION],
    SubgridDims[MAX_DIMENSION];

  /* make sure it is 3D */
  
  if (MetaData.TopGridRank != 3) {
    ENZO_VFAIL("Cannot do GalaxySimulation in %"ISYM" dimension(s)\n", MetaData.TopGridRank)
  }

  /* set default parameters */

  float GalaxySimulationGasMass,
    GalaxySimulationGalaxyMass,
    GalaxySimulationDiskTemperature,
    GalaxySimulationAngularMomentum[MAX_DIMENSION],
    GalaxySimulationUniformVelocity[MAX_DIMENSION],
    GalaxySimulationUniformDensity,
    GalaxySimulationUniformEnergy,
    GalaxySimulationAmbientDensity;

  FLOAT GalaxySimulationDiskRadius,
    GalaxySimulationDiskPosition[MAX_DIMENSION],
    GalaxySimulationDiskScaleHeightz,
    GalaxySimulationDiskScaleHeightR,
	GalaxySimulationSubgridLeft[MAX_DIMENSION],
    GalaxySimulationSubgridRight[MAX_DIMENSION];

  float GalaxySimulationInitialTemperature,
    GalaxySimulationDarkMatterConcentrationParameter,
    GalaxySimulationInflowTime,
    GalaxySimulationInflowDensity;

  int   GalaxySimulationRefineAtStart,
    GalaxySimulationUseMetallicityField;

  int GalaxySimulationB0FieldFlag;
  float GalaxySimulationB0FieldStrength;
  
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  float ZeroBField[3] = {0.0, 0.0, 0.0};

  /* Default Values */

  GalaxySimulationRefineAtStart      = FALSE;
  GalaxySimulationUseMetallicityField  = FALSE;
  GalaxySimulationInitialTemperature = 1000.0;
  GalaxySimulationDiskRadius         = 0.2;      // [Mpc]
  GalaxySimulationDiskTemperature    = 1.e4;     // [K]
  GalaxySimulationDiskScaleHeightz   = 325e-6;
  GalaxySimulationDiskScaleHeightR   = 3500e-6;
  GalaxySimulationGasMass            = 4.0e10;
  GalaxySimulationDiskTemperature    = 1000.0;
  GalaxySimulationAmbientDensity     = 1.674e-30; // [g cm-3]
  GalaxySimulationInflowTime         = -1;
  GalaxySimulationInflowDensity      = 0;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    GalaxySimulationDiskPosition[dim] = 0.5*(DomainLeftEdge[dim] +
					     DomainRightEdge[dim]);
    GalaxySimulationAngularMomentum[dim] = 0;
    GalaxySimulationUniformVelocity[dim] = 0;
  }
  GalaxySimulationUniformDensity = 1.0;
  GalaxySimulationUniformEnergy = 1.0;
  GalaxySimulationSubgridLeft[0] = GalaxySimulationSubgridLeft[1] =
  GalaxySimulationSubgridLeft[2] = 0.0;    // start of subgrid(s)
  GalaxySimulationSubgridRight[0] = GalaxySimulationSubgridRight[1] =
  GalaxySimulationSubgridRight[2] = 0.0;    // end of subgrid(s)

    
  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    
    ret = 0;
   
    ret += sscanf(line, "GalaxySimulationRefineAtStart = %"ISYM,
		  &GalaxySimulationRefineAtStart);
    ret += sscanf(line, "GalaxySimulationSubgridLeft = %"FSYM" %"FSYM" %"FSYM,
                    GalaxySimulationSubgridLeft,GalaxySimulationSubgridLeft+1,GalaxySimulationSubgridLeft+2);
    ret += sscanf(line, "GalaxySimulationSubgridRight = %"FSYM" %"FSYM" %"FSYM,
                    GalaxySimulationSubgridRight,GalaxySimulationSubgridRight+1,GalaxySimulationSubgridRight+2);
    ret += sscanf(line, "GalaxySimulationUseMetallicityField = %"ISYM,
		  &GalaxySimulationUseMetallicityField);
    ret += sscanf(line, "GalaxySimulationInitialTemperature = %"FSYM,
		  &GalaxySimulationInitialTemperature);
    ret += sscanf(line, "GalaxySimulationUniformVelocity = %"FSYM" %"FSYM" %"FSYM,
                  &GalaxySimulationUniformVelocity[0], &GalaxySimulationUniformVelocity[1],
                  &GalaxySimulationUniformVelocity[2]);
    ret += sscanf(line, "GalaxySimulationDiskRadius = %"PSYM,
		  &GalaxySimulationDiskRadius);
    ret += sscanf(line, "GalaxySimulationGasMass = %"FSYM,
		  &GalaxySimulationGasMass);
    ret += sscanf(line, "GalaxySimulationDiskPosition = %"PSYM" %"PSYM" %"PSYM, 
		  &GalaxySimulationDiskPosition[0],
		  &GalaxySimulationDiskPosition[1],
		  &GalaxySimulationDiskPosition[2]);
    ret += sscanf(line, "GalaxySimulationDiskScaleHeightz = %"PSYM,
		  &GalaxySimulationDiskScaleHeightz);
    ret += sscanf(line, "GalaxySimulationDiskScaleHeightR = %"PSYM,
		  &GalaxySimulationDiskScaleHeightR);
    ret += sscanf(line, "GalaxySimulationDiskTemperature = %"FSYM,
		  &GalaxySimulationDiskTemperature);
    ret += sscanf(line, "GalaxySimulationAmbientDensity = %"FSYM,
                    &GalaxySimulationAmbientDensity);
    ret += sscanf(line, "GalaxySimulationInflowTime = %"FSYM,
		  &GalaxySimulationInflowTime);
    ret += sscanf(line, "GalaxySimulationInflowDensity = %"FSYM,
		  &GalaxySimulationInflowDensity);
    ret += sscanf(line, "GalaxySimulationAngularMomentum = %"FSYM" %"FSYM" %"FSYM,
		  &GalaxySimulationAngularMomentum[0],
		  &GalaxySimulationAngularMomentum[1],
		  &GalaxySimulationAngularMomentum[2]);
    ret += sscanf(line, "GalaxySimulationB0FieldFlag = %"ISYM,
		  &GalaxySimulationB0FieldFlag);
    ret += sscanf(line, "GalaxySimulationB0FieldStrength = %"FSYM,
		  &GalaxySimulationB0FieldStrength);
    
    /* if the line is suspicious, issue a warning */
    
    if (ret == 0 && strstr(line, "=") && strstr(line, "GalaxySimulation") 
	&& line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file

 /* set up field names and units */

 int count = 0,mm;
 DataLabel[count++] = DensName;
 DataLabel[count++] = TEName;
 if (DualEnergyFormalism)
   DataLabel[count++] = GEName;
 DataLabel[count++] = Vel1Name;
 if(MetaData.TopGridRank > 1)
   DataLabel[count++] = Vel2Name;
 if(MetaData.TopGridRank > 2)
   DataLabel[count++] = Vel3Name;
 DataLabel[count++] = G0Name;
 DataLabel[count++] = ZetaName;
 printf("HydroMethod:%d\n",HydroMethod);
 if( HydroMethod == MHD_RK || HydroMethod == MHD_Li ){
   DataLabel[count++] = B1Name;
   DataLabel[count++] = B2Name;
   DataLabel[count++] = B3Name;
 }
 if( HydroMethod == MHD_RK){
   DataLabel[count++] = PhiName;
 }
 if (GalaxySimulationUseMetallicityField)
   DataLabel[count++] = MetalName;
 if (StarMakerTypeIaSNe)
   DataLabel[count++] = MetalIaName;
 for (i = 0; i < count; i++)
   DataUnits[i] = NULL;
 for(mm=0;mm<count;mm++)
   printf("%s\n",DataLabel[mm]);
    
if (UseMHDCT){
   MHDcLabel[0] = "Bx";
   MHDcLabel[1] = "By";
   MHDcLabel[2] = "Bz";
   MHDLabel[0] = "BxF";
   MHDLabel[1] = "ByF";
   MHDLabel[2] = "BzF";
   MHDeLabel[0] = "Ex";
   MHDeLabel[1] = "Ey";
   MHDeLabel[2] = "Ez";
   
   MHDUnits[0] = "None";
   MHDUnits[1] = "None";
   MHDUnits[2] = "None";
   MHDeUnits[0] = "None";
   MHDeUnits[1] = "None";
   MHDeUnits[2] = "None";
}

  /* set up grid */

  if (TopGrid.GridData->MyGalaxyInitializeGrid(GalaxySimulationDiskRadius,
						       ExternalGravityConstant, 
						       GalaxySimulationGasMass,
						       GalaxySimulationDiskPosition, 
						       GalaxySimulationDiskScaleHeightz,
						       GalaxySimulationDiskScaleHeightR, 
						       GalaxySimulationDiskTemperature, 
						       GalaxySimulationInitialTemperature,
						       GalaxySimulationAngularMomentum,
						       GalaxySimulationUniformVelocity,
                               GalaxySimulationAmbientDensity,
						       GalaxySimulationUseMetallicityField,
						       GalaxySimulationInflowTime,
						       GalaxySimulationInflowDensity,
							   GalaxySimulationB0FieldFlag,
							   GalaxySimulationB0FieldStrength,0)
	      == FAIL) {
      ENZO_FAIL("Error in GalaxySimulationInitialize[Sub]Grid.");
  }// end subgrid if

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (int dim = 0; dim < MetaData.TopGridRank; dim++)
      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
	float(MetaData.TopGridDims[dim]);
  }

  /* If requested, refine the grid to the desired level. */

  if (GalaxySimulationRefineAtStart) {
      /* Create as many subgrids as there are refinement levels
       needed to resolve the initial explosion region upon the start-up. */
      
      HierarchyEntry ** Subgrid;
      if (MaximumRefinementLevel > 0)
      Subgrid   = new HierarchyEntry*[MaximumRefinementLevel];
      
      /* Create new HierarchyEntries. */
      
      int lev;
      for (lev = 0; lev < MaximumRefinementLevel; lev++)
      Subgrid[lev] = new HierarchyEntry;
      
      for (lev = 0; lev < MaximumRefinementLevel; lev++) {
          
          for (dim = 0; dim < MetaData.TopGridRank; dim++)
          NumberOfSubgridZones[dim] =
          nint((GalaxySimulationSubgridRight[dim] - GalaxySimulationSubgridLeft[dim])/
               ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
                float(MetaData.TopGridDims[dim])))
          *int(POW(RefineBy, lev + 1));
          
          if (debug)
          printf("RotatingCylinder:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1,
                 NumberOfSubgridZones[0]);
          
          if (NumberOfSubgridZones[0] > 0) {
              
              /* fill them out */
              
              if (lev == 0)
              TopGrid.NextGridNextLevel  = Subgrid[0];
              Subgrid[lev]->NextGridThisLevel = NULL;
              if (lev == MaximumRefinementLevel-1)
              Subgrid[lev]->NextGridNextLevel = NULL;
              else
              Subgrid[lev]->NextGridNextLevel = Subgrid[lev+1];
              if (lev == 0)
              Subgrid[lev]->ParentGrid        = &TopGrid;
              else
              Subgrid[lev]->ParentGrid        = Subgrid[lev-1];
              
              /* compute the dimensions and left/right edges for the subgrid */
              
              for (dim = 0; dim < MetaData.TopGridRank; dim++) {
                  SubgridDims[dim] = NumberOfSubgridZones[dim] + 2*NumberOfGhostZones;
                  LeftEdge[dim]    = GalaxySimulationSubgridLeft[dim];
                  RightEdge[dim]   = GalaxySimulationSubgridRight[dim];
              }
              
              /* create a new subgrid and initialize it */
              
              Subgrid[lev]->GridData = new grid;
              Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
              Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
                                                  LeftEdge, RightEdge, 0);
              if (Subgrid[lev]->GridData->MyGalaxyInitializeGrid(GalaxySimulationDiskRadius,
                                                                 ExternalGravityConstant,
                                                                 GalaxySimulationGasMass,
                                                                 GalaxySimulationDiskPosition,
                                                                 GalaxySimulationDiskScaleHeightz,
                                                                 GalaxySimulationDiskScaleHeightR,
                                                                 GalaxySimulationDiskTemperature,
                                                                 GalaxySimulationInitialTemperature,
                                                                 GalaxySimulationAngularMomentum,
                                                                 GalaxySimulationUniformVelocity,
                                                                 GalaxySimulationAmbientDensity,
                                                                 GalaxySimulationUseMetallicityField,
                                                                 GalaxySimulationInflowTime,
                                                                 GalaxySimulationInflowDensity,
																 GalaxySimulationB0FieldFlag,
																 GalaxySimulationB0FieldStrength,0)
                  == FAIL) {
                  ENZO_FAIL("Error in InitializeUniformGrid (subgrid).");
              }
              
              /* set up the initial explosion area on the finest resolution subgrid */
              
              if (lev == MaximumRefinementLevel - 1)
              if (Subgrid[lev]->GridData->MyGalaxyInitializeGrid(GalaxySimulationDiskRadius,
                                                                 ExternalGravityConstant,
                                                                 GalaxySimulationGasMass,
                                                                 GalaxySimulationDiskPosition,
                                                                 GalaxySimulationDiskScaleHeightz,
                                                                 GalaxySimulationDiskScaleHeightR,
                                                                 GalaxySimulationDiskTemperature,
                                                                 GalaxySimulationInitialTemperature,
                                                                 GalaxySimulationAngularMomentum,
                                                                 GalaxySimulationUniformVelocity,
                                                                 GalaxySimulationAmbientDensity,
                                                                 GalaxySimulationUseMetallicityField,
                                                                 GalaxySimulationInflowTime,
                                                                 GalaxySimulationInflowDensity,
																 GalaxySimulationB0FieldFlag,
																 GalaxySimulationB0FieldStrength,0)
                  == FAIL) {
                  ENZO_FAIL("Error in RotatingCylinderInitialize[Sub]Grid.");
              }
              
          }
          else{
              printf("RotatingCylinder: single grid start-up.\n");
          }
      }
      
      
      /* set up subgrids from level 1 to max refinement level -1 */
      
      for (lev = MaximumRefinementLevel - 1; lev > 0; lev--)
      if (Subgrid[lev]->GridData->ProjectSolutionToParentGrid(
                                                              *(Subgrid[lev-1]->GridData))
          == FAIL) {
          ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
      }
      
      /* set up the root grid */
      
      if (MaximumRefinementLevel > 0) {
          if (Subgrid[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData))
              == FAIL) {
              ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
          }
      }
      else
      if (TopGrid.GridData->MyGalaxyInitializeGrid(GalaxySimulationDiskRadius,
                                                   ExternalGravityConstant,
                                                   GalaxySimulationGasMass,
                                                   GalaxySimulationDiskPosition,
                                                   GalaxySimulationDiskScaleHeightz,
                                                   GalaxySimulationDiskScaleHeightR,
                                                   GalaxySimulationDiskTemperature,
                                                   GalaxySimulationInitialTemperature,
                                                   GalaxySimulationAngularMomentum,
                                                   GalaxySimulationUniformVelocity,
                                                   GalaxySimulationAmbientDensity,
                                                   GalaxySimulationUseMetallicityField,
                                                   GalaxySimulationInflowTime,
                                                   GalaxySimulationInflowDensity,
												   GalaxySimulationB0FieldFlag,
												   GalaxySimulationB0FieldStrength,0)
          == FAIL) {
          ENZO_FAIL("Error in RotatingCylinderInitializeGrid.");
      }

  } // end: if (GalaxySimulationRefineAtStart)


 /* Write parameters to parameter output file */

 if (MyProcessorNumber == ROOT_PROCESSOR) {

   fprintf(Outfptr, "GalaxySimulationRefineAtStart      = %"ISYM"\n",
	   GalaxySimulationRefineAtStart);
   fprintf(Outfptr, "GalaxySimulationUseMetallicityField          = %"ISYM"\n",
	   GalaxySimulationUseMetallicityField);
   fprintf(Outfptr, "GalaxySimulationInitialTemperature = %"GOUTSYM"\n",
	   GalaxySimulationInitialTemperature);
   fprintf(Outfptr, "GalaxySimulationUniformVelocity    = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"\n",
	   GalaxySimulationUniformVelocity[0], GalaxySimulationUniformVelocity[1],
	   GalaxySimulationUniformVelocity[2]);
   fprintf(Outfptr, "GalaxySimulationDiskRadius = %"GOUTSYM"\n",
	   GalaxySimulationDiskRadius);
   fprintf(Outfptr, "DiskRotatingVelocity = %"GOUTSYM"\n",
	   ExternalGravityConstant);
   fprintf(Outfptr, "GalaxySimulationGasMass = %"GOUTSYM"\n",
	   GalaxySimulationGasMass);
   fprintf(Outfptr, "GalaxySimulationDiskScaleHeightz = %"GOUTSYM"\n",
	   GalaxySimulationDiskScaleHeightz);
   fprintf(Outfptr, "GalaxySimulationDiskScaleHeightR = %"GOUTSYM"\n",
	   GalaxySimulationDiskScaleHeightR);
   fprintf(Outfptr, "GalaxySimulationDiskTemperature = %"GOUTSYM"\n",
	   GalaxySimulationDiskTemperature);
   fprintf(Outfptr, "GalaxySimulationAmbientDensity = %"GOUTSYM"\n",
	   GalaxySimulationAmbientDensity);
   fprintf(Outfptr, "GalaxySimulationInflowTime = %"GOUTSYM"\n",
	   GalaxySimulationInflowTime);
   fprintf(Outfptr, "GalaxySimulationInflowDensity = %"GOUTSYM"\n",
	   GalaxySimulationInflowDensity);
   fprintf(Outfptr, "GalaxySimulationDiskPosition = ");
   WriteListOfFloats(Outfptr, MetaData.TopGridRank, GalaxySimulationDiskPosition);
   fprintf(Outfptr, "GalaxySimulationAngularMomentum = ");
   WriteListOfFloats(Outfptr, MetaData.TopGridRank, GalaxySimulationAngularMomentum);
   fprintf(Outfptr, "GalaxySimulationB0FieldFlag = %"ISYM"\n",
	   GalaxySimulationB0FieldFlag);
   fprintf(Outfptr, "GalaxySimulationB0FieldStrength = %"GOUTSYM"\n",
	   GalaxySimulationB0FieldStrength);
 }

#ifdef USE_MPI

 // BWO: this forces the synchronization of the various point source gravity
 // parameters between processors.  If this is not done, things go to pieces!

 MPI_Barrier(MPI_COMM_WORLD);
 MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
 MPI_Bcast(&PointSourceGravityConstant,1,DataType,ROOT_PROCESSOR, MPI_COMM_WORLD);
 MPI_Bcast(&PointSourceGravityCoreRadius,1,DataType,ROOT_PROCESSOR, MPI_COMM_WORLD);

#endif

 return SUCCESS;

}
