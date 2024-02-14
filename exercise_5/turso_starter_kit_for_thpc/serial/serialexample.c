#include <stdio.h>
#include <stdlib.h>

int main(int argc,char *argv[])
{

  printf("host: %s\n",getenv("HOSTNAME"));
  exit(0);

}
