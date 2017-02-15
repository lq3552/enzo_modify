#include <stdio.h>
void auto_show_flags(FILE *fp) {
   fprintf (fp,"\n");
   fprintf (fp,"CPP = /apps/compilers/gcc/5.2.0/bin/cpp\n");
   fprintf (fp,"CC  = /apps/mpi/gcc/5.2.0/openmpi/1.10.2/bin/mpicc\n");
   fprintf (fp,"CXX = /apps/mpi/gcc/5.2.0/openmpi/1.10.2/bin/mpic++\n");
   fprintf (fp,"FC  = /apps/compilers/gcc/5.2.0/bin/gfortran\n");
   fprintf (fp,"F90 = /apps/compilers/gcc/5.2.0/bin/gfortran\n");
   fprintf (fp,"LD  = /apps/mpi/gcc/5.2.0/openmpi/1.10.2/bin/mpic++\n");
   fprintf (fp,"\n");
   fprintf (fp,"DEFINES = -DLINUX -DH5_USE_16_API   -D__max_subgrids=100000 -D__max_baryons=30 -D__max_cpu_per_node=8 -D__memory_pool_size=100000 -DINITS64 -DLARGE_INTS -DCONFIG_PINT_8 -DIO_32    -DUSE_MPI   -DCONFIG_PFLOAT_8 -DCONFIG_BFLOAT_8  -DUSE_HDF5_GROUPS   -DTRANSFER   -DNEW_GRID_IO -DFAST_SIB      -DENZO_PERFORMANCE  -DUSE_GRACKLE  -DSAB\n");
   fprintf (fp,"\n");
   fprintf (fp,"INCLUDES = -I/apps/hdf5/1.8.9/include  -I/apps/gcc/5.2.0/libz/1.2.8/include        -I/ufrc/tan/pg3552/local/include    -I.\n");
   fprintf (fp,"\n");
   fprintf (fp,"CPPFLAGS = -P -traditional \n");
   fprintf (fp,"CFLAGS   =  -O2\n");
   fprintf (fp,"CXXFLAGS =  -O2\n");
   fprintf (fp,"FFLAGS   = -fno-second-underscore -ffixed-line-length-132 -O2\n");
   fprintf (fp,"F90FLAGS = -fno-second-underscore -O2\n");
   fprintf (fp,"LDFLAGS  =  -O2\n");
   fprintf (fp,"\n");
   fprintf (fp,"LIBS     = -L/apps/hdf5/1.8.9/lib -lhdf5 -lz  -lgfortran  -L/apps/gcc/5.2.0/libz/1.2.8/lib -lz        -L/ufrc/tan/pg3552/local/lib -lgrackle\n");
   fprintf (fp,"\n");
}
