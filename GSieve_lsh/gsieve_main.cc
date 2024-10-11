#include <unistd.h>
#include <getopt.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include <NTL/LLL.h>
#include "gsieve.h"
#include "math.h"

NTL_CLIENT

const char* USAGE =
  "Usage: %s [OPTION]\n"
  "-f FILE_NAME:\tA file containing the input basis. "
    "If not set, read from stdin\n"
  "-g GOAL_NORM:\tIf the algorithm finds a vector v, "
    "with |v|^2 = GOAL_NORM it will stop.\n"
    "\t\tEither way the algorithm terminates after a lot of "
    "collisions are found.\n"
  "-v:\t\tVerbose mode, print additional information.\n"
  "-h:\t\tPrint this help.\n";

int main(int argc, char** argv) {
  char* input_file_name = NULL;
  long goal_norm = 0;
  bool flag_verbose = false;
  int option;
  while((option = getopt (argc, argv, "f:g:hv")) != -1) {
    switch (option) {
      case 'v':
        flag_verbose = true;
        break;
      case 'f':
        input_file_name = optarg;
        break;
      case 'g':
        goal_norm = atol(optarg);
        if (goal_norm < 0)
          goal_norm = 0;
        break;
      case 'h':
      case '?':
        fprintf (stderr, USAGE, argv[0]);
        return -1;
    }
  }
  mat_ZZ B;
  ifstream input_file(input_file_name);
  if (input_file.is_open()) {
    input_file >> B;
    input_file.close();
  } else {
    cin >> B;
  }
  // Reduce the Basis, BKZ with window 20
  G_BKZ_FP(B, 0.99, 20);
  KleinSampler sampler;
  GSieve gsieve;
  gsieve.Init(B, &sampler);
  gsieve.SetGoalSqrNorm(goal_norm);
  gsieve.SetVerbose(flag_verbose);
  gsieve.Start();
}

