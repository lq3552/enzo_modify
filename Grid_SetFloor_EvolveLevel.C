/***********************************************************************
/
/  GRID CLASS (SET INTERNAL ENERGY FLOOR)
/
/  written by: Peng Wang
/  date:       October, 2007
/  modified1:
/
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "TopGridData.h"
#include "Grid.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::SetFloor_EvolveLevel()
{

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }
  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1, TimeUnits, 
    VelocityUnits, CriticalDensity = 1, BoxLength = 1, MagneticUnits;
  double MassUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, 1.0);// 1.0 -> time

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, 
    B1Num, B2Num, B3Num, HMNum, H2INum, H2IINum;
  int ibx,iby,ibz;
  
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) { //Here G0Num is not universally used so...
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }

	iden = DensNum;
	ietot =TENum;
    if(DualEnergyFormalism){
		ieint = GENum;
	}
	ivx = Vel1Num;
	ivy = Vel2Num;
	ivz = Vel3Num;
	if (UseMHDCT){
		ibx = B1Num;
		iby = B2Num;
		ibz = B3Num;
	}
#if 1
  FILE *fp;
  float vx, vy, vz, v2, eint,eint0,etot0, emin,eminJ, rho, bx,by,bz,b2;
  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {

        int igrid = (k * GridDimension[1] + j) * GridDimension[0] + i;
        rho = BaryonField[iden][igrid];
        vx = BaryonField[ivx][igrid];
        vy = BaryonField[ivy][igrid];
        vz = BaryonField[ivz][igrid];
        v2 = vx*vx + vy*vy + vz*vz;
		if (UseMHDCT){
			bx = BaryonField[ibx][igrid];
			by = BaryonField[iby][igrid];
			bz = BaryonField[ibz][igrid];
			b2 = bx*bx + by*by + bz*bz;
		}
        
		// Internal energy limiter
		if(HydroMethod == Zeus_Hydro){
          eint = BaryonField[ietot][igrid]; // I'll try but this may be wrong
        }
		else if (DualEnergyFormalism){
          eint = BaryonField[ieint][igrid];
        }
		else if (UseMHDCT){
		  eint = BaryonField[ietot][igrid] - 0.5*v2 - 0.5*b2/rho;
		}
        else {
          eint = BaryonField[ietot][igrid] - 0.5*v2;
        }
        
        if(UseGasTemperatureFloor){ // Fixed Gas Temperature Floor
            emin = GasTemperatureFloor/TemperatureUnits/((Gamma-1.0)*Mu); // Mu should be adapted for Multispecies in the future
		}
        else{
            emin = 0.0;
        }
        /* 4xgrid Jeans length support! */
        eminJ = 0.0;//4.0*0.48999*rho*pow(CellWidth[0][0],2)/(Gamma*(Gamma-1.0));
		eint0 = eint;
		etot0 = BaryonField[ietot][igrid];
        eint = max(eint, emin);
        eint = max(eint, eminJ);
		if(HydroMethod == Zeus_Hydro)
			BaryonField[ietot][igrid] = eint;
		else if (UseMHDCT)
			BaryonField[ietot][igrid] = eint + 0.5*v2 + 0.5*b2/rho;
		else 
			BaryonField[ietot][igrid] = eint + 0.5*v2;
        if (DualEnergyFormalism) {
          BaryonField[ieint][igrid] = eint;
        }
		if (etot0 < 0 ||eint0 < 0|| rho < 0){ // error message
			fp = fopen("trace_outlier.txt","a");
			fprintf(fp,"%d\t%d\t%d\t%f\t%f\t%f\t%d\t%d\t\n",i,j,k,eint0*((Gamma-1.0)*Mu)*TemperatureUnits,etot0,rho,Time,ieint,ietot);
		}
		
		// alven speed limiter, limit rho... I feel bad ...
		if (HydroMethod == MHD_Li) {
			const float ca_min = MaximumAlvenSpeed;
			float ca = sqrt(b2)/sqrt(rho);
			if (ca > ca_min) {
				BaryonField[ietot][igrid] -= 0.5*b2/rho;
				float rho1 = b2/pow(ca_min,2);
				BaryonField[iden][igrid] = rho1;
				BaryonField[ietot][igrid] += 0.5*b2/rho1;
				printf("density floor set based on MaximumAlvenSpeed: (%"GSYM" %"GSYM" %"GSYM"), rho: %"GSYM"->%"GSYM"\n", CellLeftEdge[0][i],
					CellLeftEdge[1][j], CellLeftEdge[2][k], rho*DensityUnits, rho1*DensityUnits);
			}
		}
      } // for i
    } // for j
  } // for k
#endif
    

  
  return SUCCESS;
}
