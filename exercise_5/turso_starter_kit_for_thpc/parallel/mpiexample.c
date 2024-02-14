#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>

int main(int argc,char *argv[])
{
  const int tag=50;
  int id,ntasks,source_id,dest_id,rc,i;
  MPI_Status status;
  int msg[2],pnlen;
  char pname[MPI_MAX_PROCESSOR_NAME];

  rc=MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    printf("MPI initialization failed\n");
    exit(1);
  }
  rc=MPI_Comm_size(MPI_COMM_WORLD,&ntasks);
  rc=MPI_Comm_rank(MPI_COMM_WORLD,&id);
  rc=MPI_Get_processor_name(pname,&pnlen);

  if (id != 0) {
    msg[0]=id;
    msg[1]=ntasks;
    dest_id=0;
    rc=MPI_Send(msg,2,MPI_INT,dest_id,tag,MPI_COMM_WORLD);
  } else {
    for (i=1; i < ntasks; i++) {
      rc=MPI_Recv(msg,2,MPI_INT,MPI_ANY_SOURCE,tag,MPI_COMM_WORLD,&status);
      source_id=status.MPI_SOURCE;
      printf("message: %d %d sender: %d\n",msg[0],msg[1],source_id);
    }
  }

  printf("\nhost %s\n",pname);
  rc=MPI_Finalize();
  exit(0);
}
