#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

double evaluate_f_of_x(double x) {
        int k;
        float coeffs[100];
        double f_of_x;
        for (size_t i = 0; i < k; i++) {
          f_of_x = f_of_x + (coeffs[i] + pow(x, i));
        }
        return f_of_x;
      }

double integrate(double a, double b, int m) {
         double x;
         double f_of_x;
         double h = ( b - a ) / m;
         double integral = 0.0;
         for (int i = 0; i <= m; i++) {
           x = a + i * h;
           f_of_x = evaluate_f_of_x(x);
           if ( ( i == 0 ) || ( i == m ) )
             integral += 0.5 * f_of_x;
           else
             integral += f_of_x;
         }
         integral *= h;
         return integral;
       }

int isNumberCharacter (char asciiCharacter) {
  if ((asciiCharacter >= 48 && asciiCharacter <= 57) ||
    asciiCharacter == 45 || asciiCharacter == 46) {
    return 1;
  } else {
    return 0;
  }
}

int isEndLineCharacter (char asciiCharacter) {
  if (asciiCharacter == 10 || asciiCharacter == 13) {
    return 1;
  } else {
    return 0;
  }
}

int getSingleIntFromLine(FILE *file) {
  char cache[10];
  int size = 0;
  char readAsciiCharacter = 0;

  while ((readAsciiCharacter = getc(file)) != EOF && isNumberCharacter(readAsciiCharacter)) {
    cache[size] = readAsciiCharacter;
    size++;
  }
  return atoi(cache);
}

void cleanTable (char *table) {
  for (int i = 0; i < 100; ++i) {
    table[i] = 0;
  }
}

int getTableOfArguments(FILE *file, float *array) {

  char cache[100];
  int arraySize = 0,
      cacheSize = 0,
      asciiCharacter = 0;

  while ((asciiCharacter = getc(file)) != EOF || asciiCharacter == 10) {
          if (isNumberCharacter(asciiCharacter)) {
            cache[cacheSize] = asciiCharacter;
            cacheSize ++;
          }
          if (asciiCharacter == 32 || isEndLineCharacter(asciiCharacter)) {
            array[arraySize] = atof(cache);
            arraySize ++;
            cacheSize = 0;
            cleanTable(cache);
            if (isEndLineCharacter(asciiCharacter)) {
              break;
            }
         }
       }
       return arraySize;
}

void getInterval(FILE *file, float *start, float *end) {

  char cache[100];
  float startInterval = 0.0,
        endInterval = 0.0;
  int secondArgument = 0,
      cacheSize = 0,
      asciiCharacter = 0;
  while ((asciiCharacter = getc(file)) != EOF || asciiCharacter == 10) {
          if (isNumberCharacter(asciiCharacter)) {
            cache[cacheSize] = asciiCharacter;
            cacheSize ++;
          }
          if (asciiCharacter == 32 || isEndLineCharacter(asciiCharacter)) {
            if (secondArgument) {
              endInterval = atof(cache);
              break;
            } else {
              secondArgument = 1;
              startInterval = atof(cache);
              cacheSize = 0;
              cleanTable(cache);
            }
         }
       }
  *start = startInterval;
  *end = endInterval;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    char parametersTypeCache;
    char *parametersSingleValueCache;
    float coeffs[100];
    float startInterval = 0,
          endInterval = 0;
    int asciiCharacter = 0,
        coeffsSize = 0;
    int k, integration;

    FILE *file;
    file = fopen("test.txt", "r");
    if (file) {
        while ((asciiCharacter = getc(file)) != EOF){
          if ((asciiCharacter >= 'a' && asciiCharacter <= 'z') || asciiCharacter == 32) {
            if (asciiCharacter == 32) {
              switch(parametersTypeCache) {
                case 'e':
                  printf("degree:\n");
                  k = getSingleIntFromLine(file);
                  printf("%i\n", k);
                break;
                case 's':
                  printf("coeffs:\n");
                  coeffsSize = getTableOfArguments(file, coeffs);
                  for (size_t i = 0; i < coeffsSize; i ++) {
                      printf("%lf\n", coeffs[i]);
                  }
                break;
                case 'l':
                  printf("interval:\n");
                  getInterval(file, &startInterval, &endInterval);
                  printf("%f\n", startInterval);
                  printf("%f\n", endInterval);
                break;
                case 'n':
                  printf("integration:\n");
                  integration = getSingleIntFromLine(file);
                  printf("%i\n", integration);
                break;
                default:
                  printf("none");
              }
            }
            parametersTypeCache = asciiCharacter;
         }
        }
    }
    //float calka;
    //calka = integrate(startInterval, endInterval, integration);
    //printf ("%f\n", calka);

    int pack_size;
    char *pack_buff;
    int pack_position;
    int all_floats = k + 3;
    int all_ints = 2;

    if (rank == 0) {
      int tmp_pack_size;

      MPI_Pack_size(all_ints, MPI_INT, MPI_COMM_WORLD, &tmp_pack_size);
      pack_size = tmp_pack_size;

      MPI_Pack_size(all_floats, MPI_FLOAT, MPI_COMM_WORLD, &tmp_pack_size);
      pack_size += tmp_pack_size;

      pack_buff = new char [pack_size];

      pack_position = 0;
      MPI_Pack(&k, all_ints, MPI_INT, pack_buff, pack_size, &pack_position, MPI_COMM_WORLD);
      MPI_Pack(&coeffs[0], all_floats, MPI_FLOAT, pack_buff, pack_size, &pack_position, MPI_COMM_WORLD);

      MPI_Bcast(%pack_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(pack_buff, pack_size, MPI_PACKED, 0, MPI_COMM_WORLD);
    } else {
      MPI_Bcast(&pack_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
      pack_buff = new char [pack_size];
      MPI_Bcast(pack_buff, pack_size, MPI_PACKED, 0, MPI_COMM_WORLD);

      pack_position = 0;
      MPI_Unpack(pack_buff, pack_size, &pack_position, &k, all_ints, MPI_INT, MPI_COMM_WORLD);
      MPI_Unpack(pack_buff, pack_size, &pack_position, &coeffs[0], all_floats, MPI_FLOAT, MPI_COMM_WORLD);
    }

    //dekompozycja
    int p = m / n;
    int r = m % n;
    int my_number_of_subintervals;
    int my_first_midpoint;
    if (rank < r) {
      my_number_of_subintervals = p + 1;
      my_first_midpoint = rank * (p + 1);
    } else {
      my_number_of_subintervals = p;
      my_first_midpoint = r * (p + 1) + (rank - r) * p;
    }

    MPI_Finalize();
    return 0;
}
