/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A GALAXY SIMULATION)
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:  Elizabeth Tasker, Feb, 2004
/  modified1:  Elizabeth Tasker, Oct, 2006 (tidied up)
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "hydro_rk/EOS.h"

#define Mpc (3.0856e24)         //Mpc [cm] 
#define SolarMass (1.989e33)    //Solar Mass [g]
#define GravConst (6.67e-8)     //Gravitational Constant [cm3g-1s-2]
#define pi (3.14159)
#define mh (1.674e-24)           //Mass of Hydrogen [g]
#define kboltz (1.381e-16)      //Boltzmann's Constant [ergK-1]

//!!! That is important, notice the units!!
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

/* Internal routines */
float mygauss_mass(FLOAT r, FLOAT z, FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT inv [3][3], float DiskDensity, FLOAT ScaleHeightR, FLOAT ScaleHeightz, FLOAT cellwidth);
void myrot_to_disk(FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT &xrot, FLOAT &yrot, FLOAT &zrot, FLOAT inv [3][3]);

static float DensityUnits, LengthUnits, TemperatureUnits = 1, TimeUnits, VelocityUnits;

int grid::MyGalaxyInitializeGrid(FLOAT DiskRadius,
					 float CircularVelocity,
					 float GasMass,
					 FLOAT DiskPosition[MAX_DIMENSION], 
					 FLOAT ScaleHeightz,
					 FLOAT ScaleHeightR, 
					 float DiskTemperature,
					 float InitialTemperature,
					 float AngularMomentum[MAX_DIMENSION],
					 float UniformVelocity[MAX_DIMENSION],
                     float AmbientDensity,
					 int UseMetallicityField, 
					 float GalaxySimulationInflowTime,
					 float GalaxySimulationInflowDensity,
					 int B0Flag,
					 float B0Strength,
					 int level)
{
 /* declarations */
	
  int dim, i, j, k, m, field, disk, size=1, MetalNum, MetalIaNum, vel;
 int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
   DINum, DIINum, HDINum, B1Num, B2Num, B3Num, PhiNum,G0Num,ZetaNum;
 float DiskDensity, DiskVelocityMag;
 float Scale[3];

  /* Record field type for FindField */

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism)
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  vel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1) 
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;
  FieldType[G0Num = NumberOfBaryonFields++] = G0;
  FieldType[ZetaNum = NumberOfBaryonFields++] = Zeta;

  if(HydroMethod == MHD_RK || HydroMethod == MHD_Li){
	  UseMHD = 1;
	  FieldType[B1Num=NumberOfBaryonFields++] = Bfield1;
	  FieldType[B2Num=NumberOfBaryonFields++] = Bfield2;
	  FieldType[B3Num=NumberOfBaryonFields++] = Bfield3;
  }
  if(HydroMethod == MHD_RK){
	  FieldType[PhiNum=NumberOfBaryonFields++] = PhiField;
  }

  if (MultiSpecies) { // Include Chemistry
    FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
    if (MultiSpecies > 1) {
      FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
      FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
      FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
    }
    if (MultiSpecies > 2) {
      FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
    }
  }

  if (UseMetallicityField) // Include Metallicity-2
    FieldType[MetalNum = NumberOfBaryonFields++] = Metallicity; /* fake it with metals */
  if (StarMakerTypeIaSNe)
    FieldType[MetalIaNum = NumberOfBaryonFields++] = MetalSNIaDensity;

 /* Return if this doesn't concern us. */

 if (ProcessorNumber != MyProcessorNumber) 
   return SUCCESS;

 /* Set various units. */

 float CriticalDensity = 1, BoxLength = 1, mu = Mu;//mu0=0.6 (herewith an important question!!
 FLOAT a, dadt, ExpansionFactor = 1;
 
 /* Set up units */
 if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
	ENZO_FAIL("Error in GetUnits.\n");
 }
 float MagneticUnits = sqrt(DensityUnits*4.0*M_PI)*VelocityUnits; 

 /* Set up inflow */
 if (GalaxySimulationInflowTime > 0.0){
   TimeActionType[0] = 2;
   TimeActionParameter[0] = GalaxySimulationInflowDensity*DensityUnits;
   TimeActionTime[0] = GalaxySimulationInflowTime*1e9/TimeUnits;
 }

 /* compute size of fields */

 for (dim = 0; dim < GridRank; dim++){
    size *= GridDimension[dim];
	Scale[dim]=(GridRightEdge[dim]-GridLeftEdge[dim])/(GridDimension[dim]-2*NumberOfGhostZones);
 }
 /* allocate fields */
 
 this->AllocateGrids();

/*for (field = 0; field < NumberOfBaryonFields; field++)
   if (BaryonField[field] == NULL)
     BaryonField[field] = new float[size];*/

 int ijk,dd = NumberOfGhostZones;
 float density, density1, density2, dens1, Velocity[MAX_DIMENSION],FUV_G0,CR_Zeta,
   temperature, temp1,phi_G,*phi_S,phi_T;
 phi_S = new float[size];
 FLOAT r, x, y = 0, z = 0;
 float q = 0.7;
 int n = 0, mm = 0;
 float sig1 = 9.0e+5, sig2 = 6.0e+5; //[cms-1]
 FLOAT Q1 = 1.5, Q2 = 20.;
 FLOAT RJ = 4.1*Mpc/1000., R0 = 8.5*Mpc/1000.; // scale R of J_FUV; zero point [cm]
 float g0 = 1.7; 
 float kappa;
 FILE *testfile,*SG;

 /* load self-gravity field */
 if (!(SG = fopen("SG.in","r")))
     ENZO_FAIL("Fail to load self-gravity field!\n");
 printf("NumberofGhostZone = %d\nGridDimension = %d %d %d\n",dd,GridDimension[0],GridDimension[1],GridDimension[2]);
 for (i = dd; i < GridDimension[0]-dd; i++)
     for (j = dd; j < GridDimension[1]-dd; j++)
         for (k = dd; k < GridDimension[2]-dd; k++){
             ijk = i*GridDimension[1]*GridDimension[2]+j*GridDimension[2]+k;
             fscanf(SG,"%lf\n",&phi_S[ijk]);
         } //end load self-gravity field
 fclose(SG);
    
 density1 = sig1/(2*M_PI*GravConst*Q1*ScaleHeightz*Mpc); //[cgs]
 density2 = sig2/(2*M_PI*GravConst*Q2*ScaleHeightz*Mpc); //[cgs]
 //EOSSoundSpeed=sqrt(InitialTemperature/TemperatureUnits/mu); //code unit
 printf("soundspeed=%lf\n",EOSSoundSpeed);
 printf("ExternalGravityConstant=%lf\n",ExternalGravityConstant);
 
 printf("UseMHDCT = %d\n",UseMHDCT);
 /* MHD specific field initialize */
 if(UseMHDCT){
	 if (B0Flag == 1){ // Uniform field
		float UniformField[3] = {B0Strength,0,0}; ///!!! Will become another parameter in later version!
		for ( int field=0; field < 3; field++ ){
			for (i=0; i<MagneticSize[field]; i++ ){
				MagneticField[field][i] = UniformField[field]/MagneticUnits;
			}
		}
	 } //end if 1
	 else if (B0Flag == 2){ //Toroidal field
		float X,Y,Z,R,rc = 0.002*Mpc/LengthUnits;
		int index,field;
		for (field=0;field<3;field++){
			for( k=0;k<ElectricDims[field][2];k++){
				for( j=0;j<ElectricDims[field][1];j++){
					for( i=0;i<ElectricDims[field][0];i++){
						index = i+ElectricDims[field][0]*(j+ElectricDims[field][1]*k);
						if (field==2){
							X = (i-GridStartIndex[0])*Scale[0]-0.5-0.5*Scale[0];
							Y = (j-GridStartIndex[1])*Scale[1]-0.5-0.5*Scale[1];
							Z = (k-GridStartIndex[2])*Scale[2]-0.5-0.5*Scale[2];
							R = sqrt(POW(X,2)+POW(Y,2));
							if (R > rc){
								ElectricField[field][index] = -B0Strength*R/MagneticUnits/POW(cosh(Z*LengthUnits/ScaleHeightz/Mpc),2);
							}
							else{
								ElectricField[field][index] = (-2*B0Strength/rc*R*R
                                                               +B0Strength/(rc*rc)*R*R*R)/MagneticUnits/POW(cosh(Z*LengthUnits/ScaleHeightz/Mpc),2);
							}
						}
						else
							ElectricField[field][index] = 0.0;
					}
				}
			}
		}
		if( this->MHD_Curl(GridStartIndex, GridEndIndex, 0) == FAIL ){
			fprintf(stderr," error occored in MHD_Curl\n"); 
			return FAIL;
		}
	 } // end if 2
     else if (B0Flag == 3){ //Toroidal field with ambient treatment
         float X,Y,Z,R,rc = 0.002*Mpc/LengthUnits,rh1 = 0.010*Mpc/LengthUnits,rh2 = 0.012*Mpc/LengthUnits,A,B,C,D;
         int index,field;
         A = (0.1*(0. - 1.0*B0Strength*POW(rh1,9)*POW(rh2,2) +
                                   37.0*B0Strength*POW(rh1,8)*POW(rh2,3) + 1.*B0Strength*POW(rh1,11)*POW(rh2,3) -
                                   25.0*B0Strength*POW(rh1,7)*POW(rh2,4) - 2.0*B0Strength*POW(rh1,10)*POW(rh2,4) -
                                   47.0*B0Strength*POW(rh1,6)*POW(rh2,5) + 10.0*B0Strength*POW(rh1,9)*POW(rh2,5) +
                                   26.0*B0Strength*POW(rh1,5)*POW(rh2,6) - 18.0*B0Strength*POW(rh1,8)*POW(rh2,6) +
                                   10.0*B0Strength*POW(rh1,4)*POW(rh2,7) - 1.0*B0Strength*POW(rh1,7)*POW(rh2,7) +
                                   20.0*B0Strength*POW(rh1,6)*POW(rh2,8) - 10.0*B0Strength*POW(rh1,5)*POW(rh2,9)))/
         (-1.*POW(rh1,9)*POW(rh2,2) + 1.*POW(rh1,8)*POW(rh2,3) + 1.*POW(rh1,11)*POW(rh2,3) + 2.*POW(rh1,7)*POW(rh2,4) -
          2.*POW(rh1,10)*POW(rh2,4) - 2.*POW(rh1,6)*POW(rh2,5) + 1.*POW(rh1,9)*POW(rh2,5) - 1.*POW(rh1,5)*POW(rh2,6) +
          1.*POW(rh1,4)*POW(rh2,7) - 1.*POW(rh1,7)*POW(rh2,7) + 2.*POW(rh1,6)*POW(rh2,8) - 1.*POW(rh1,5)*POW(rh2,9));
         B = (-1.*(1.8*B0Strength*POW(rh1,6)*POW(rh2,3) -
                   5.4*B0Strength*POW(rh1,5)*POW(rh2,4) + 1.8*B0Strength*POW(rh1,8)*POW(rh2,4) +
                   10.8*B0Strength*POW(rh1,3)*POW(rh2,6) -
                   10.8*B0Strength*POW(rh1,6)*POW(rh2,6) + 1.8*B0Strength*POW(rh1,9)*POW(rh2,6) -
                   5.4*B0Strength*POW(rh1,2)*POW(rh2,7) + 10.8*B0Strength*POW(rh1,5)*POW(rh2,7) -
                   3.6*B0Strength*POW(rh1,8)*POW(rh2,7) - 5.4*B0Strength*rh1*POW(rh2,8) +
                   5.4*B0Strength*POW(rh1,4)*POW(rh2,8) +
                   3.6*B0Strength*POW(rh2,9) - 10.8*B0Strength*POW(rh1,3)*POW(rh2,9) +
                   3.6*B0Strength*POW(rh1,6)*POW(rh2,9) + 3.6*B0Strength*POW(rh1,2)*POW(rh2,10) -
                   1.8*B0Strength*POW(rh1,5)*POW(rh2,10)))/
         ((1.*POW(rh1,2)*rh2 - 3.*rh1*POW(rh2,2) + 1.*POW(rh1,4)*POW(rh2,2) + 2.*POW(rh2,3) - 1.*POW(rh1,3)*POW(rh2,3))*
          (-1.*POW(rh1,5)*POW(rh2,2) + 1.*POW(rh1,4)*POW(rh2,3) + 1.*POW(rh1,7)*POW(rh2,3) + 2.*POW(rh1,3)*POW(rh2,4) -
           2.*POW(rh1,6)*POW(rh2,4) - 2.*POW(rh1,2)*POW(rh2,5) + 1.*POW(rh1,5)*POW(rh2,5) - 1.*rh1*POW(rh2,6) + 1.*POW(rh2,7) -
           1.*POW(rh1,3)*POW(rh2,7) + 2.*POW(rh1,2)*POW(rh2,8) - 1.*rh1*POW(rh2,9)));
         C = (-1.*(-0.9*B0Strength*POW(rh1,9)*POW(rh2,4) +
                   3.6*B0Strength*POW(rh1,8)*POW(rh2,5) - 1.8*B0Strength*POW(rh1,11)*POW(rh2,5)
                    - 1.8*B0Strength*POW(rh1,7)*POW(rh2,6) +
                   4.5*B0Strength*POW(rh1,10)*POW(rh2,6) - 0.9*B0Strength*POW(rh1,13)*POW(rh2,6) -
                   9.0*B0Strength*POW(rh1,6)*POW(rh2,7) + 3.6*B0Strength*POW(rh1,9)*POW(rh2,7) +
                   0.9*B0Strength*POW(rh1,12)*POW(rh2,7) + 10.8*B0Strength*POW(rh1,5)*POW(rh2,8) -
                   16.2*B0Strength*POW(rh1,8)*POW(rh2,8) + 3.6*B0Strength*POW(rh1,11)*POW(rh2,8) +
                   5.4*B0Strength*POW(rh1,4)*POW(rh2,9) + 1.8*B0Strength*POW(rh1,7)*POW(rh2,9) -
                   3.6*B0Strength*POW(rh1,10)*POW(rh2,9) - 12.6*B0Strength*POW(rh1,3)*POW(rh2,10) +
                   21.6*B0Strength*POW(rh1,6)*POW(rh2,10) - 5.4*B0Strength*POW(rh1,9)*POW(rh2,10) +
                   1.8*B0Strength*POW(rh1,2)*POW(rh2,11) - 9.0*B0Strength*POW(rh1,5)*POW(rh2,11) +
                   5.4*B0Strength*POW(rh1,8)*POW(rh2,11) + 4.5*B0Strength*rh1*POW(rh2,12) -
                   12.6*B0Strength*POW(rh1,4)*POW(rh2,12) + 3.6*B0Strength*POW(rh1,7)*POW(rh2,12) -
                   1.8*B0Strength*POW(rh2,13) + 7.2*B0Strength*POW(rh1,3)*POW(rh2,13) -
                   3.6*B0Strength*POW(rh1,6)*POW(rh2,13) + 2.7*B0Strength*POW(rh1,2)*POW(rh2,14) -
                   0.9*B0Strength*POW(rh1,5)*POW(rh2,14) - 1.8*B0Strength*rh1*POW(rh2,15) + 0.9*B0Strength*POW(rh1,4)*POW(rh2,15)))/
         ((1.*POW(rh1,2) - 1.*POW(rh2,2))*(1.*POW(rh1,2)*rh2 - 1.*POW(rh2,3))*
          (1.*POW(rh1,2)*rh2 - 3.*rh1*POW(rh2,2) + 1.*POW(rh1,4)*POW(rh2,2) + 2.*POW(rh2,3) - 1.*POW(rh1,3)*POW(rh2,3))*
          (-1.*POW(rh1,5)*POW(rh2,2) + 1.*POW(rh1,4)*POW(rh2,3) + 1.*POW(rh1,7)*POW(rh2,3) + 2.*POW(rh1,3)*POW(rh2,4) -
           2.*POW(rh1,6)*POW(rh2,4) - 2.*POW(rh1,2)*POW(rh2,5) + 1.*POW(rh1,5)*POW(rh2,5) - 1.*rh1*POW(rh2,6) + 1.*POW(rh2,7) -
           1.*POW(rh1,3)*POW(rh2,7) + 2.*POW(rh1,2)*POW(rh2,8) - 1.*rh1*POW(rh2,9)));
         D = (-1.*(0. + 0.45*B0Strength*POW(rh1,10)*POW(rh2,5) - 1.8*B0Strength*POW(rh1,9)*POW(rh2,6) +
                   0.45*B0Strength*POW(rh1,12)*POW(rh2,6) + 0.9*B0Strength*POW(rh1,8)*POW(rh2,7) -
                   0.9*B0Strength*POW(rh1,11)*POW(rh2,7) + 4.5*B0Strength*POW(rh1,7)*POW(rh2,8) -
                   0.9*B0Strength*POW(rh1,10)*POW(rh2,8) - 5.4*B0Strength*POW(rh1,6)*POW(rh2,9) +
                   2.7*B0Strength*POW(rh1,9)*POW(rh2,9) - 2.7*B0Strength*POW(rh1,5)*POW(rh2,10) +
                   6.3*B0Strength*POW(rh1,4)*POW(rh2,11) - 2.7*B0Strength*POW(rh1,7)*POW(rh2,11) -
                   0.9*B0Strength*POW(rh1,3)*POW(rh2,12) + 0.9*B0Strength*POW(rh1,6)*POW(rh2,12) -
                   2.25*B0Strength*POW(rh1,2)*POW(rh2,13) + 0.9*B0Strength*POW(rh1,5)*POW(rh2,13) +
                   0.9*B0Strength*rh1*POW(rh2,14) - 0.45*B0Strength*POW(rh1,4)*POW(rh2,14)))/
         ((1.*POW(rh1,2) - 1.*POW(rh2,2))*(1.*POW(rh1,2)*rh2 - 1.*POW(rh2,3))*
          (1.*POW(rh1,2)*rh2 - 3.*rh1*POW(rh2,2) + 1.*POW(rh1,4)*POW(rh2,2) + 2.*POW(rh2,3) - 1.*POW(rh1,3)*POW(rh2,3))*
          (-1.*POW(rh1,5)*POW(rh2,2) + 1.*POW(rh1,4)*POW(rh2,3) + 1.*POW(rh1,7)*POW(rh2,3) + 2.*POW(rh1,3)*POW(rh2,4) -
           2.*POW(rh1,6)*POW(rh2,4) - 2.*POW(rh1,2)*POW(rh2,5) + 1.*POW(rh1,5)*POW(rh2,5) - 1.*rh1*POW(rh2,6) + 1.*POW(rh2,7) -
           1.*POW(rh1,3)*POW(rh2,7) + 2.*POW(rh1,2)*POW(rh2,8) - 1.*rh1*POW(rh2,9)));
         for (field=0;field<3;field++){
             for( k=0;k<ElectricDims[field][2];k++){
                 for( j=0;j<ElectricDims[field][1];j++){
                     for( i=0;i<ElectricDims[field][0];i++){
                         index = i+ElectricDims[field][0]*(j+ElectricDims[field][1]*k);
                         if (field==2){
                             X = (i-GridStartIndex[0])*Scale[0]-0.5-0.5*Scale[0];
                             Y = (j-GridStartIndex[1])*Scale[1]-0.5-0.5*Scale[1];
                             Z = (k-GridStartIndex[2])*Scale[2]-0.5-0.5*Scale[2];
                             R = sqrt(POW(X,2)+POW(Y,2));
                             if (R > rc && R < rh1)
                                 ElectricField[field][index] = -B0Strength*R/MagneticUnits/POW(cosh(Z*LengthUnits/ScaleHeightz/Mpc),2);
                             else if (R <= rc)
                                 ElectricField[field][index] = (-2*B0Strength/rc*R*R
                                                                +B0Strength/(rc*rc)*R*R*R)/MagneticUnits/POW(cosh(Z*LengthUnits/ScaleHeightz/Mpc),2);
                             else if (R >= rh1 && R < rh2)
								 ElectricField[field][index] = -(A*R+B*R*R+C*R*R*R+D*R*R*R*R)/MagneticUnits/POW(cosh(Z*LengthUnits/ScaleHeightz/Mpc),2);
							 else
								 ElectricField[field][index] = -B0Strength*R/10./MagneticUnits/POW(cosh(Z*LengthUnits/ScaleHeightz/Mpc),2);
                         }
                         else
                             ElectricField[field][index] = 0.0;
                     }
                 }
             }
         }
         if( this->MHD_Curl(GridStartIndex, GridEndIndex, 0) == FAIL ){
             fprintf(stderr," error occored in MHD_Curl\n"); 
             return FAIL;
         }
     } // endif 3
	 this->CenterMagneticField();
	 float *DivB=NULL;
	 MHD_Diagnose("Post Initialize Grid", DivB);
 } // end if(UseMHDCT)

 /* loop begin grids */
 for (k = 0; k < GridDimension[2]; k++)
   for (j = 0; j < GridDimension[1]; j++)
     for (i = 0; i < GridDimension[0]; i++, n++) {
    ijk = (i-dd)*GridDimension[1]*GridDimension[2]+(j-dd)*GridDimension[2]+(k-dd);
	/* Compute position */
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	if (GridRank > 1)
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	if (GridRank > 2)
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

	density = AmbientDensity/DensityUnits; //1 or less? Here I use the fixed unit
	temperature = temp1 = InitialTemperature;
         
	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  Velocity[dim] = 0;

	/* Find distance from center. */

	r = sqrt(POW(fabs(x-DiskPosition[0]), 2) +
		 POW(fabs(y-DiskPosition[1]), 2));
	r = max(r, 0.1*CellWidth[0][0]);

	if (r < DiskRadius) {

	  FLOAT xpos, ypos, zpos, zheight, drad;
	  float CellMass;
	  FLOAT xhat[3];
	  FLOAT yhat[3];



	  /* Loop over dims if using Zeus (since vel's face-centered). */

	  for (dim = 0; dim < 1+(HydroMethod == Zeus_Hydro ? GridRank : 0);
	       dim++) { //dim 0 --- sys motion; dim 1-3: surface velocity for viscosity calculation!!!!

	    /* Compute position. */

	    xpos = x-DiskPosition[0] - 
	      (dim == 1 ? 0.5*CellWidth[0][0] : 0.0);
	    ypos = y-DiskPosition[1] -
	      (dim == 2 ? 0.5*CellWidth[1][0] : 0.0);
	    zpos = z-DiskPosition[2] -
	      (dim == 3 ? 0.5*CellWidth[2][0] : 0.0);
	    
	    /* Compute z and r_perp (AngularMomentum is angular momentum 
	       and must have unit length). We don't need to care about that
		   z-axis as default rotating axis !!! */    

	    /* magnitude of z = r.L in L direction */

	    zheight = AngularMomentum[0]*xpos + 
	              AngularMomentum[1]*ypos +
	              AngularMomentum[2]*zpos;

	    /* position in plane of disk */

	    xhat[0] = xpos - zheight*AngularMomentum[0];
	    xhat[1] = ypos - zheight*AngularMomentum[1];
	    xhat[2] = zpos - zheight*AngularMomentum[2]; //Just 2D, no azimu-D
	    drad = sqrt(xhat[0]*xhat[0] + xhat[1]*xhat[1] + xhat[2]*xhat[2]);


	    /* Normalize the vector r_perp = unit vector pointing along plane of disk */

	    xhat[0] = xhat[0]/drad;
	    xhat[1] = xhat[1]/drad;
	    xhat[2] = xhat[2]/drad;

	    /* Find another vector perpendicular to r_perp and AngularMomentum */

	    yhat[0] = AngularMomentum[1]*xhat[2] - AngularMomentum[2]*xhat[1];
	    yhat[1] = AngularMomentum[2]*xhat[0] - AngularMomentum[0]*xhat[2];
	    yhat[2] = AngularMomentum[0]*xhat[1] - AngularMomentum[1]*xhat[0];

	    /* generate rotation matrix */
	    FLOAT inv[3][3],temp;
	    int i,j;
	    
	    // matrix of basis vectors in coordinate system defined by the galaxy
	    inv[0][0] = xhat[0]; inv[0][1] = yhat[0]; inv[0][2] = AngularMomentum[0];
	    inv[1][0] = xhat[1]; inv[1][1] = yhat[1]; inv[1][2] = AngularMomentum[1];
	    inv[2][0] = xhat[2]; inv[2][1] = yhat[2]; inv[2][2] = AngularMomentum[2];
	    
	    // Matrix is orthogonal by construction so inverse = transpose (??Needed?)
	    for (i=0;i<3;i++)
	      for (j=i+1;j<3;j++)
		{
		  temp = inv[i][j];
		  inv[i][j] = inv[j][i];
		  inv[j][i] = temp;
		}

        /* Compute velocity magnitude. This assumes Logarithm Potential*/
          //Diskdensity=
	    DiskVelocityMag = CircularVelocity*r/sqrt(POW(drad,2)+POW(ExternalGravityRadius,2)); //Code units
          //(r>ExternalGravityRadius?CircularVelocity*drad/sqrt(POW(r,2)+POW(ExternalGravityRadius,2)):0.0)
        kappa = sqrt(2)*DiskVelocityMag/r*sqrt(1+r/DiskVelocityMag*(CircularVelocity*POW(ExternalGravityRadius,2)/POW(POW(ExternalGravityRadius,2)+POW(r,2),3.0/2.0)))/TimeUnits; //check the units! Important
	    if (dim == 0)
	      { // So the density is an "average" of a whole cell ...
		//CellMass = mygauss_mass(drad*LengthUnits,zheight*LengthUnits, xpos*LengthUnits, ypos*LengthUnits, zpos*LengthUnits, inv,
				     // DiskDensity*DensityUnits,ScaleHeightR*Mpc, ScaleHeightz*Mpc, CellWidth[0][0]*LengthUnits);
        // I think here, POW(CellWidth[0][0],3) is agressive, then grid must be cubic
		//dens1 = CellMass/POW(CellWidth[0][0]*LengthUnits,3)/DensityUnits;
          dens1=((r*LengthUnits>0.002*Mpc && r*LengthUnits<0.01*Mpc) ? kappa*density1/POW(cosh((z-DiskPosition[2])*LengthUnits/ScaleHeightz/Mpc),2) : kappa*density2/POW(cosh((z-DiskPosition[2])*LengthUnits/ScaleHeightz/Mpc),2))/DensityUnits;
	      }
        

	    if (dens1 < density){
	      break;
		}

	    /* Compute velocty: L x r_perp. */

	    if (dim == 0 || dim == 1)
	      Velocity[0] = DiskVelocityMag*(AngularMomentum[1]*xhat[2] -
					     AngularMomentum[2]*xhat[1]);
	    if (dim == 0 || dim == 2)
	      Velocity[1] = DiskVelocityMag*(AngularMomentum[2]*xhat[0] -
					     AngularMomentum[0]*xhat[2]);
	    if (dim == 0 || dim == 3)
	      Velocity[2] = DiskVelocityMag*(AngularMomentum[0]*xhat[1] -
					     AngularMomentum[1]*xhat[0]);
	    
	  } // end: loop over dims

	   	    
	    /* If the density is larger than the background (or the previous
	       disk), then set as gas cell. */
	  if (dens1 > density) {
	    density = dens1;
		FUV_G0 = (r*LengthUnits<4.0*Mpc/1000. ? g0*exp(-(4.0*Mpc/1000.-R0)/RJ) : g0*exp(-(r*LengthUnits-R0)/RJ));
	    CR_Zeta = (r*LengthUnits<10.0*Mpc/1000.? (5e-16-(5e-16-5e-17)/(0.01*Mpc)*r) : 5e-17);
		if (temp1 == InitialTemperature)
	      temp1 = DiskTemperature;
	    temperature = temp1;
	  }
      else {
          //set density field to balance background gravity
          phi_G = 0.5*pow(CircularVelocity*VelocityUnits,2)
              *log(pow(ExternalGravityRadius*LengthUnits,2)*pow(ExternalGravityRadius*LengthUnits,-2)
              +pow(r*LengthUnits,2)*pow(ExternalGravityRadius*LengthUnits,-2)
              +pow((z-DiskPosition[2])*LengthUnits,2)*pow(ExternalGravityRadius*LengthUnits,-2)/pow(q,2))
              /pow(LengthUnits/TimeUnits,2);// Gravitational potential
          if (i>=dd&&j>=dd&&k>=dd&&i<GridDimension[0]-dd&&j<GridDimension[1]-dd&&k<GridDimension[2]-dd)
              phi_T = phi_G-phi_S[ijk]/pow(LengthUnits/TimeUnits,2);
          else
              phi_T = phi_G; // here!!!
          FUV_G0 = 0.01; //fake one
		  CR_Zeta = 5e-17;//fake?? Should be zero...!!!..!!!...!!
		  density = AmbientDensity/DensityUnits*exp(-(mu*TemperatureUnits*phi_T/InitialTemperature));
      }
	} // end: if (r < DiskRadius)
    else{
        //set density field to balance background gravity
        phi_G = 0.5*pow(CircularVelocity*VelocityUnits,2)
            *log(pow(ExternalGravityRadius*LengthUnits,2)*pow(ExternalGravityRadius*LengthUnits,-2)
            +pow(r*LengthUnits,2)*pow(ExternalGravityRadius*LengthUnits,-2)
            +pow((z-DiskPosition[2])*LengthUnits,2)*pow(ExternalGravityRadius*LengthUnits,-2)/pow(q,2))
            /pow(LengthUnits/TimeUnits,2);// Gravitational potential
          if (i>=dd&&j>=dd&&k>=dd&&i<GridDimension[0]-dd&&j<GridDimension[1]-dd&&k<GridDimension[2]-dd)
              phi_T = phi_G-phi_S[ijk]/pow(LengthUnits/TimeUnits,2);
          else
              phi_T = phi_G; // here!!!
		  FUV_G0 = 0.01;
		  CR_Zeta = 5e-17;
		  density = AmbientDensity/DensityUnits*exp(-(mu*TemperatureUnits*phi_T/InitialTemperature));
    }
	/* Set density. */

	BaryonField[0][n] = density;
	BaryonField[G0Num][n] = FUV_G0;	
	BaryonField[ZetaNum][n] = CR_Zeta;	
	/* Fake metallicity. */
	if (UseMetallicityField)
	    BaryonField[MetalNum][n] = 0.0;
	if (StarMakerTypeIaSNe)
	    BaryonField[MetalIaNum][n] = 1.0e-10;
	
	/* Fake B fields. */
	if(HydroMethod == MHD_RK){
	    BaryonField[B1Num][n] = 0.0;
	    BaryonField[B2Num][n] = 0.0;
	    BaryonField[B3Num][n] = 0.1;
        BaryonField[PhiNum][n] = 0.0;
	}

	/* Set Velocities. */

	for (dim = 0; dim < GridRank; dim++)
	  BaryonField[vel+dim][n] = Velocity[dim] + UniformVelocity[dim];

	/* Set energy (thermal and then total if necessary). */

	BaryonField[1][n] = temperature/TemperatureUnits/
                           ((Gamma-1.0)*mu);

	if (DualEnergyFormalism)
	  BaryonField[2][n] = BaryonField[1][n]; 
	
	if (HydroMethod != Zeus_Hydro)
	  for (dim = 0; dim < GridRank; dim++)
	    BaryonField[1][n] += 0.5*POW(BaryonField[vel+dim][n], 2);
    if (UseMHD)
        BaryonField[1][n] += 0.5*(POW(BaryonField[B1Num][n],2)+POW(BaryonField[B2Num][n],2)+POW(BaryonField[B3Num][n],2))/BaryonField[0][n];
    if (BaryonField[1][n]<0)
	  printf("n = %d  T = %g   e = %g\n", n, temperature,BaryonField[1][n]);
   } // end loop over grid (Baryonfield initialize
      
 return SUCCESS;

}


// Computes the total mass in a given cell by integrating the density profile using 5-point Gaussian quadrature (I can change it if using Tasker&Tan 2009)
float mygauss_mass(FLOAT r, FLOAT z, FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT inv [3][3], float DiskDensity, FLOAT ScaleHeightR, FLOAT ScaleHeightz, FLOAT cellwidth)
{
  
  FLOAT EvaluationPoints [5] = {-0.90617985,-0.53846931,0.0,0.53846931,0.90617985};
  FLOAT Weights [5] = {0.23692689,0.47862867,0.56888889,0.47862867,0.23692689};
  FLOAT xResult [5];
  FLOAT yResult [5];
  float Mass = 0;
  FLOAT xrot,yrot,zrot;
  int i,j,k;

  for (i=0;i<5;i++)
    {
      xResult[i] = 0.0;
      for (j=0;j<5;j++)
	{
	  yResult[j] = 0.0;
	  for (k=0;k<5;k++)
	    {
	      myrot_to_disk(xpos+EvaluationPoints[i]*cellwidth/2.0,ypos+EvaluationPoints[j]*cellwidth/2.0,zpos+EvaluationPoints[k]*cellwidth/2.0,xrot,yrot,zrot,inv);
	      yResult[j] += cellwidth/2.0*Weights[k]*PEXP(-sqrt(POW(xrot,2)+POW(yrot,2))/ScaleHeightR)/POW(cosh(zrot/(2.0*ScaleHeightz)),2);
	    }
	  xResult[i] += cellwidth/2.0*Weights[j]*yResult[j];
	}
      Mass += cellwidth/2.0*Weights[i]*xResult[i];
    }  
  Mass *= DiskDensity;
  return Mass;
}

//Finds coordinates in rotated coordinate system
void myrot_to_disk(FLOAT xpos, FLOAT ypos, FLOAT zpos, FLOAT &xrot, FLOAT &yrot, FLOAT &zrot, FLOAT inv [3][3])
{
  xrot = xpos*inv[0][0] + ypos*inv[0][1] + zpos*inv[0][2];
  yrot = xpos*inv[1][0] + ypos*inv[1][1] + zpos*inv[1][2];
  zrot = xpos*inv[2][0] + ypos*inv[2][1] + zpos*inv[2][2];
}
