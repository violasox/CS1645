#include <stdio.h>
#include <mpi.h>


int main (int argc, char* argv[])
{
  int my_rank,len,np,i;
  char pname[256];


  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
  MPI_Comm_size(MPI_COMM_WORLD,&np);
  MPI_Get_processor_name(pname,&len);

  printf("PROCESS %d of %d: Running on %s\n",my_rank+1,np,pname);


  MPI_Finalize();
}


