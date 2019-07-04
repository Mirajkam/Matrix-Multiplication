
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <vector>

using namespace std;

unsigned long *matrixMult(unsigned long *m1, unsigned long *m2, const int M,
                          const int N, const int K);

unsigned long *matrixMult(unsigned long *m1, unsigned long *m2, const int M,
                          const int N, const int K) {

  unsigned long *result = new unsigned long[M * N]; // r1 * c2

  unsigned long *temp = new unsigned long[N * K]; // c2 * r2

  // M: row of result matrix and row of matrix 1 - 1
  // N: col of result matrix and col of matrix 2 - 2
  // K: col of matrix 1 and row of matirx 2 - 1
#pragma omp parallel for collapse(2)
  for (int j = 0; j < N; ++j) {
    for (int k = 0; k < K; ++k) {
      // temp[j][k] = m2[k][j];
      temp[j * K + k] = m2[k * N + j];
    }
  }

#pragma omp parallel for collapse(2)
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      result[i * N + j] = 0;
      for (int k = 0; k < K; ++k) {
        // result[i][j] += m1[i][k]*temp[j][k];
        result[i * N + j] += m1[i * K + k] * m2[j * K + k];
      }
    }
  }

  int acc00 = 0;
  int acc01 = 0;
  int acc02 = 0;
  int acc03 = 0;

  for (int i = 0; i < M; i += 2) {
    for (int j = 0; j < N; j += 2) {
      result[i * N + j] = 0;
      acc00 = acc01 = acc02 = acc03 = 0;

      for (int k = 0; k < K; ++k) {
        // result[i][j] += m1[i][k]*temp[j][k];

        result[i * N + j] += m1[i * K + k] * m2[j * K + k];
      }
    }
  }

  return result;
}

int main(int argc, char *argv[]) {
  int R1 = atoi(argv[1]);
  int C1 = atoi(argv[2]);
  int R2 = atoi(argv[3]);
  int C2 = atoi(argv[4]);

  unsigned long *m1 = new unsigned long[R1 * C1];
  unsigned long *m2 = new unsigned long[R2 * C2];

  cout << "Matrix 1:\n";
  for (int i = 0; i < R1; i++) {
    for (int j = 0; j < C1; j++) {
      m1[i * C1 + j] = rand() % 1000 + 1;
    }
    // cout << "\n";
  }

  cout << "Matrix 2:\n";
  for (int i = 0; i < R2; i++) {
    for (int j = 0; j < C2; j++) {
      m2[i * C2 + j] = rand() % 1000 + 1;
    }
    // cout << "\n";
  }

  unsigned long *result = matrixMult(m1, m2, R1, C2, R2);
  /* cout << "\n     Matrix Result:\n";
   for (int i = 0; i < R1; i++) {
     for (int j = 0; j < C2; j++) {
       cout << (result[i * C2 + j]) << " ";
     }
     cout << "\n";
   } */

  return 0;
}
