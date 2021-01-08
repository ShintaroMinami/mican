#include "mican.h"

void print_progress(char *progress_name, int numerator, int denominator) {

  float progress;
  progress = (float)(100 * numerator / denominator);
  fprintf(stderr, "%18s ... %4d /%4d (%3.0f%%) in       sec.\n", progress_name,
          numerator, denominator, progress);
  fprintf(stderr, "\033[1A");

  return;
}
