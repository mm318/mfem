#ifndef CMFEM_EXAMPLES_COMMON_H
#define CMFEM_EXAMPLES_COMMON_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int cmfem_parse_string_option(int argc, char *argv[], int *index,
                                     const char *short_name,
                                     const char *long_name,
                                     const char **value)
{
   if (strcmp(argv[*index], short_name) == 0 ||
       strcmp(argv[*index], long_name) == 0)
   {
      if (*index + 1 >= argc)
      {
         fprintf(stderr, "Missing value for %s\n", argv[*index]);
         return 0;
      }
      *value = argv[++(*index)];
      return 1;
   }
   return -1;
}

static int cmfem_parse_int_option(int argc, char *argv[], int *index,
                                  const char *short_name,
                                  const char *long_name,
                                  int *value)
{
   if (strcmp(argv[*index], short_name) == 0 ||
       strcmp(argv[*index], long_name) == 0)
   {
      if (*index + 1 >= argc)
      {
         fprintf(stderr, "Missing value for %s\n", argv[*index]);
         return 0;
      }
      *value = atoi(argv[++(*index)]);
      return 1;
   }
   return -1;
}

static int cmfem_parse_double_option(int argc, char *argv[], int *index,
                                     const char *short_name,
                                     const char *long_name,
                                     double *value)
{
   if (strcmp(argv[*index], short_name) == 0 ||
       strcmp(argv[*index], long_name) == 0)
   {
      if (*index + 1 >= argc)
      {
         fprintf(stderr, "Missing value for %s\n", argv[*index]);
         return 0;
      }
      *value = atof(argv[++(*index)]);
      return 1;
   }
   return -1;
}

static int cmfem_parse_bool_option(int *index, char *argv[],
                                   const char *enable_short_name,
                                   const char *enable_long_name,
                                   const char *disable_short_name,
                                   const char *disable_long_name,
                                   int *value)
{
   if (strcmp(argv[*index], enable_short_name) == 0 ||
       strcmp(argv[*index], enable_long_name) == 0)
   {
      *value = 1;
      return 1;
   }
   if (strcmp(argv[*index], disable_short_name) == 0 ||
       strcmp(argv[*index], disable_long_name) == 0)
   {
      *value = 0;
      return 1;
   }
   return -1;
}

static int cmfem_uniform_refinement_levels(double max_elements,
                                           int num_elements,
                                           int dim)
{
   double levels;
   if (num_elements <= 0 || dim <= 0 || max_elements <= (double)num_elements)
   {
      return 0;
   }
   levels = floor(log(max_elements / (double)num_elements) / log(2.0) /
                  (double)dim);
   return levels > 0.0 ? (int)levels : 0;
}

#endif
