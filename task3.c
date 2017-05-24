//
//  task3.c
//  zad3
//
//  Created by Kamil Moreński on 16.05.2017.
//  Copyright © 2017 Kamil Moreński. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

double evaluate_f_of_x(int k, float coeffs[], double x){
  double sum = 0.0;
  double powOfx;
  for (size_t i = 0; i < k + 1; i++) { //z uwzglednieneim zerowego wsploczynnika
    powOfx = pow(x, k);
    sum += powOfx * coeffs[k];
  }
  return sum;
}
double analytically(int k, float coeffs[], double x){
  double sum = 0.0;
  double powOfx;
  for (size_t i = 0; i < k + 1; i++) {
    sum += coeffs * pow(x, k + 1) / (k + 1);
  }
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
    int k, integration;
    if (rank == 0) {
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
                printf("degree: ");
                k = getSingleIntFromLine(file);
                printf("%i\n", k);
              break;
              case 's':
                printf("coeffs:\n");
                coeffsSize = getTableOfArguments(file, coeffs);
                for (size_t i = 0; i < coeffsSize; i ++) {
                    printf("%f ", coeffs[i]);
                }
              break;
              case 'l':
                printf("\ninterval: ");
                getInterval(file, &startInterval, &endInterval);
                printf("%f ", startInterval);
                printf("%f\n", endInterval);
              break;
              case 'n':
                printf("integration: ");
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
  }

        int pack_size;
        char* pack_buff;
        int pack_position;
        int all_floats = k + 3;
        int all_ints = 2;

        if (rank == 0) {
          int tmp_pack_size;

          MPI_Pack_size(all_ints, MPI_INT, MPI_COMM_WORLD, &tmp_pack_size);
          pack_size = tmp_pack_size;

          MPI_Pack_size(all_floats, MPI_FLOAT, MPI_COMM_WORLD, &tmp_pack_size);
          pack_size += tmp_pack_size;
          printf("pack_size %d\n", pack_size);

          pack_buff = malloc(pack_size * sizeof(char));

          pack_position = 0;
          MPI_Pack(&k, 1, MPI_INT, pack_buff, pack_size, &pack_position, MPI_COMM_WORLD);
          MPI_Pack(&coeffs[0], k + 1, MPI_FLOAT, pack_buff, pack_size, &pack_position, MPI_COMM_WORLD);
          MPI_Pack(&startInterval, 1, MPI_FLOAT, pack_buff, pack_size, &pack_position, MPI_COMM_WORLD);
          MPI_Pack(&endInterval, 1, MPI_FLOAT, pack_buff, pack_size, &pack_position, MPI_COMM_WORLD);
          MPI_Pack(&integration, 1, MPI_INT, pack_buff, pack_size, &pack_position, MPI_COMM_WORLD);

          MPI_Bcast(&pack_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
          MPI_Bcast(pack_buff, pack_size, MPI_PACKED, 0, MPI_COMM_WORLD);
        } else {
          MPI_Bcast(&pack_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
          //printf("rank%i pack_size %i\n", rank, pack_size);

          pack_buff = malloc(pack_size * sizeof(char));
          MPI_Bcast(pack_buff, pack_size, MPI_PACKED, 0, MPI_COMM_WORLD);

          pack_position = 0;
          MPI_Unpack(pack_buff, pack_size, &pack_position, &k, 1, MPI_INT, MPI_COMM_WORLD);
          MPI_Unpack(pack_buff, pack_size, &pack_position, &coeffs[0], k + 1, MPI_FLOAT, MPI_COMM_WORLD);
          MPI_Unpack(pack_buff, pack_size, &pack_position, &startInterval, 1, MPI_FLOAT, MPI_COMM_WORLD);
          MPI_Unpack(pack_buff, pack_size, &pack_position, &endInterval, 1, MPI_FLOAT, MPI_COMM_WORLD);
          MPI_Unpack(pack_buff, pack_size, &pack_position, &integration, 1, MPI_INT, MPI_COMM_WORLD);
          printf("Rank: %i ", rank);
          printf("Degree: %i ", k);
          printf("Coeffs: ");
          for (size_t i = 0; i < k + 1; i ++) {
              printf("%f ", coeffs[i]);
          }
          printf("startInterval %f ", startInterval);
          printf("endInterval %f ", endInterval);
          printf("integration %i \n", integration);
        }

        //dekompozycja
        int p = integration / size; //wynik dzielenia
        int r = integration % size; //reszta
        int my_number_of_subintervals; //liczba podprzedziałów
        int my_first_midpoint; //pierwszy punkt środkowy
        if (rank < r) {
          my_number_of_subintervals = p + 1;
          my_first_midpoint = rank * (p + 1);
        } else {
          my_number_of_subintervals = p;
          my_first_midpoint = r * (p + 1) + (rank - r) * p;
        }

        float height = (endInterval - startInterval) / integration;
        float lr = startInterval + my_first_midpoint * height;
        float ur = lr + my_number_of_subintervals * height;

        printf("Rank %i ", rank);
        printf("l(r) %f u(r) %f \n", lr, ur);

        //calkowanie
        double sum = 0.0;
        for (size_t i = 0; i < my_number_of_subintervals; i++) {
          double xi = lr + i * height;
          double xiplus1 = lr + (i + 1) * height;
          sum += evaluate_f_of_x(k, coeffs, xi) + evaluate_f_of_x(k, coeffs, xiplus1);
        }
        sum *= height;
        sum /= 2;
        //printf("%f\n", sum);

        double sumsum_numerically;
        double sum_analytically;
        MPI_Reduce(&sum, &sumsum_numerically, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        sum_analytically = analytically(k, coeffs, endInterval) - analytically(k, coeffs, startInterval);
        if (rank == 0) {
          printf("Wynik całkowania numerycznego%f\n", sumsum_numerically);
          printf("Wynik całkowania analitycznego%f\n", sum_analytically);
        }

    MPI_Finalize();
    return 0;
}
