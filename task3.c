#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

double evaluate_f_of_x(double alfa, double k)
      {
        double x;
        for (size_t i = 0; i < k; i++) {
          x = alfa * pow(x, k);
        }
      }

double integrate(double a, double b, int m)
       {
         double x;
         double f_of_x;
         double h = ( b - a ) / m;
         double integral = 0.0;
         for (int i = 0; i <= m; i++)
         {
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
  } else{
    return 0;
  }
}

int isEndLineCharacter (char asciiCharacter) {
  if (asciiCharacter == 10 || asciiCharacter == 13) {
    return 1;
  } else{
    return 0;
  }
}

int getSingleIntFromLine(FILE *file) {
  char cache[10];
  int size = 0;
  char readAsciiCharacter = 0;

  while ((readAsciiCharacter = getc(file)) != EOF && isNumberCharacter(readAsciiCharacter)){
    cache[size] = readAsciiCharacter;
    size++;
  }
  return atoi(cache);
}

void cleanTable (char *table) {
  for (int i = 0; i < 100; ++i)
  {
    table[i] = 0;
  }
}

int getTableOfArguments(FILE *file, float *array) {

  char cache[100];
  int arraySize = 0,
      cacheSize = 0,
      asciiCharacter = 0;

  while ((asciiCharacter = getc(file)) != EOF || asciiCharacter == 10){
          if (isNumberCharacter(asciiCharacter)){
            cache[cacheSize] = asciiCharacter;
            cacheSize ++;
          }
          if (asciiCharacter == 32 || isEndLineCharacter(asciiCharacter)){
            array[arraySize] = atof(cache);
            arraySize ++;
            cacheSize = 0;
            cleanTable(cache);
            if (isEndLineCharacter(asciiCharacter)){
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
  while ((asciiCharacter = getc(file)) != EOF || asciiCharacter == 10){
          if (isNumberCharacter(asciiCharacter)){
            cache[cacheSize] = asciiCharacter;
            cacheSize ++;
          }
          if (asciiCharacter == 32 || isEndLineCharacter(asciiCharacter)){
            if (secondArgument){
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

int main(int argc, char **argv)
{
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

    FILE *file;
    file = fopen("test.txt", "r");
    if (file) {
        while ((asciiCharacter = getc(file)) != EOF){
          if ((asciiCharacter >= 'a' && asciiCharacter <= 'z') || asciiCharacter == 32) {
            if (asciiCharacter == 32)
            {
              switch(parametersTypeCache)
              {
                case 'e':
                  printf("degree:\n");
                  printf("%i\n", getSingleIntFromLine(file));
                break;
                case 's':
                  printf("coeffs:\n");
                  coeffsSize = getTableOfArguments(file, coeffs);
                  int i;
                  for (i = 0; i < coeffsSize; i ++) {
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
                  printf("%i\n", getSingleIntFromLine(file));
                break;
                default:
                  printf("none");
              }
            }
            parametersTypeCache = asciiCharacter;
         }
        }
    }
    float calka;
    calka = integrate(startInterval, endInterval, getSingleIntFromLine(file));

    MPI_Finalize();
    return 0;
}
