#include "bcsr.h"
#include "logger.h"

#include <errno.h>
#include <fcntl.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

int main (int argc, char **argv)
{
  char *MM_filename = NULL;
  char *matlab_filename = NULL;
  bool linear_bin = false;
  bool log_bin = false;
  bool verbose = false;

  int spectral_bounds_method = 0;
  int numberBins = 10;
  double bin_min = 0;
  double bin_max = 1;

  int c;
  const char *short_options = "hw:d:blvm:";
  const struct option long_options[] = {
    { "help",         no_argument,        NULL, 'h' },
    { "write-MM",     required_argument,  NULL, 'w' },
    { "write-matlab", required_argument,  NULL, 'd' },
    { "linear-bin",   no_argument,        NULL, 'b' },
    { "log-bin",      no_argument,        NULL, 'l' },
    { "verbose",      no_argument,        NULL, 'v' },
    { "method",       required_argument,  NULL, 'm' },
    { NULL, 0, NULL, 0 }
  };

  initializeLogger();

  while((c = getopt_long(argc, argv, short_options, long_options,
          NULL)) != -1)
  {
    switch(c)
    {
      case 'h':
        printf("Usage: convertBCSR [options] FILE\n");
        printf("\n");
        printf("FILE is a matrix stored in BCSR binary format.\n");
        printf("\n");
        printf("{ -h | --help }               This help\n");
        printf("{ -w | --write-MM } FILE      Write matrix in MatrixMarket format to FILE\n");
        printf("{ -d | --write-matlab } FILE  Write matrix in matlab format to FILE\n");
        printf("{ -b | --log-bin }            Bin the matrix elements by magnited (log)\n");
        printf("{ -l | --linear-bin }         Bin the matrix elements by magnited (linear)\n");
        printf("{ -v | --verbose }            Print BCSR index lists\n");
        printf("{ -m | --method } N           Use method N for spectral bounds\n");
        exit(0);
        break;

      case 'w':
        MM_filename = strdup(optarg);
        break;

      case 'd':
        matlab_filename = strdup(optarg);
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

      case 'v':
        verbose = true;
        break;

      case 'm':
        spectral_bounds_method = strtol(optarg, NULL, 10);
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
    A.getSpectralBounds(spectral_bounds_method, &F_min, &F_max);
    printf("spectral bounds, method %d: [ %e, %e ]\n",
        spectral_bounds_method, F_min, F_max);
    A.toStr(verbose);

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

    if(matlab_filename != NULL)
    {
      int fd = open(matlab_filename, O_CREAT | O_EXCL | O_WRONLY, 00644);

      if(fd == -1)
      {
        if(errno == EEXIST)
        {
          ABORT("file \"%s\" already exists\n", matlab_filename);
        }

        else
        {
          ABORT("error accessing file: %s\n", strerror(errno));
        }
      }

      else
      {
        int M, N;
        double *ADense;
        A.toDense(&M, &N, &ADense);

        FILE *fstream = fdopen(fd, "w");

        fprintf(fstream, "%% %dx%d matrix\n", M, N);
        fprintf(fstream, "A = [\n");
        for(int i = 0; i < M; i++)
        {
          for(int j = 0; j < N; j++)
          {
            fprintf(fstream, " % e", ADense[i*N+j]);
          }
          fprintf(fstream, "\n");
        }
        fprintf(fstream, "];\n");

        if(fclose(fstream) != 0)
        {
          ABORT("error closing file\n");
        }
      }
    }

    exit(0);
  }

  else
  {
    ABORT("missing filename\n");
  }
}
