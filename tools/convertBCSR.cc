#include "bcsr.h"
#include "logger.h"

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main (int argc, char **argv)
{
  char *MM_filename = NULL;
  bool linear_bin = false;
  bool log_bin = false;

  int numberBins = 10;
  double bin_min = 0;
  double bin_max = 1;

  int c;
  const char *short_options = "hw:bl";
  const struct option long_options[] = {
    { "help",         no_argument,        NULL, 'h' },
    { "write-MM",     required_argument,  NULL, 'w' },
    { "linear-bin",   no_argument,        NULL, 'b' },
    { "log-bin",      no_argument,        NULL, 'l' },
    { NULL, 0, NULL, 0 }
  };

  initializeLogger();

  while((c = getopt_long(argc, argv, short_options, long_options,
          NULL)) != -1)
  {
    switch(c)
    {
      case 'h':
        printf("Usage:\n");
        printf("\n");
        printf("{ -h | --help }           This help\n");
        printf("{ -w | --write-MM } FILE  Write matrix in MatrixMarket format to FILE\n");
        printf("{ -b | --linear-bin }     Bin the matrix elements by magnited (linear)\n");
        printf("{ -b | --log-bin }        Bin the matrix elements by magnited (log)\n");
        exit(0);
        break;

      case 'w':
        MM_filename = strdup(optarg);
        break;

      case 'b':
        if(log_bin)
        {
          ABORT("log-bin already set\n");
        }
        linear_bin = true;
        break;

      case 'l':
        if(linear_bin)
        {
          ABORT("linear-bin already set\n");
        }
        log_bin = true;
        break;

      default:
        ABORT("unknown command line option\n");
        break;
    }
  }

  if(optind < argc)
  {
    BCSR A(argv[optind]);
    double F_min, F_max;
    A.getSpectralBounds(0, &F_min, &F_max);
    printf("spectral bounds Gershgorin: [ %e, %e ]\n", F_min, F_max);
    //A.getSpectralBounds(1, &F_min, &F_max);
    //printf("spectral bounds eigensolve: [ %e, %e ]\n", F_min, F_max);
    A.toStr();

    if(linear_bin || log_bin)
    {
      int *count = new int[numberBins+1];
      memset(count, 0, sizeof(int)*(numberBins+1));

      if(log_bin && bin_min < 1e-10)
      {
        bin_min = 1e-10;
      }

      double *bin = new double[numberBins+1];
      for(int i = 0; i < numberBins; i++)
      {
        if(log_bin)
        {
          bin[i] = log(bin_min)+i*(log(bin_max)-log(bin_min))/(double) numberBins;
        }

        else if(linear_bin)
        {
          bin[i] = bin_min+i*(bin_max-bin_min)/(double) numberBins;
        }
      }

      if(log_bin)
      {
        bin[numberBins] = log(bin_max);
      }

      else if(linear_bin)
      {
        bin[numberBins] = bin_max;
      }

      int outside = 0;
      int zeros = 0;
      for(int i = 0; i < A.getNumberNonZero(); i++)
      {
        bool binned = false;
        double Aij = fabs(A.getElement(i));

        if(Aij < 1e-20)
        {
          zeros++;
        }

        else
        {
          if(log_bin)
          {
            Aij = log(Aij);
          }

          for(int j = 0; j < numberBins; j++)
          {
            if(bin[j] <= Aij && Aij < bin[j+1])
            {
              binned = true;
              count[j]++;
              break;
            }
          }

          if(!binned)
          {
            outside++;
          }
        }
      }

      for(int i = 0; i < numberBins; i++)
      {
        if(log_bin)
        {
          printf("[ %e, %e ) = %d\n", exp(bin[i]), exp(bin[i+1]), count[i]);
        }

        else if(linear_bin)
        {
          printf("[ %e, %e ) = %d\n", bin[i], bin[i+1], count[i]);
        }
      }
      printf("%d elements fell outside the bins, counted %d zeros\n", outside, zeros);
    }

    if(MM_filename != NULL)
    {
      A.toMM(MM_filename);
    }

    exit(0);
  }

  else
  {
    ABORT("missing filename\n");
  }
}
