/**
 * @file rdtest.c
 * @brief Test reading the tables of scale radius and velocity created by 
 * the rdtables program.
 * @author Tomas E. Tecce
 */
#include "allvars_test.h"
#include "proto_test.h"

int main(int argc, char **argv)
{
  double chalo, mb, md, zcurr, lambda;
  double rscale, vc;


  if (argc < 2) {
    fprintf(stderr, "Error (main): no path to table files given - Exit\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
  }
  load_tables_header(*(argv+1));

  printf("Redshift?\n");
  scanf("%le", &zcurr);
  printf("Halo concentration?\n");
  scanf("%le", &chalo);
  printf("Disc mass fraction?\n");
  scanf("%le", &md);
  printf("Spin parameter?\n");
  scanf("%le", &lambda);
  printf("Bulge mass fraction?\n");
  scanf("%le", &mb);
  printf("\n");

  load_rscale_tables(zcurr, *(argv+1), *(argv+1));

  rscale = get_data_from_table(Table0_rscale, Table1_rscale, zcurr, chalo, 
			       mb, md, lambda);

  vc     = get_data_from_table(VcTable0_rscale, VcTable1_rscale, zcurr, 
			       chalo, mb, md, lambda);

  printf("Scale radius: %g * (M_200/rho_crit)^1/3\n", rscale);
  printf("Circular velocity: %g * V_200\n", vc);

  return EXIT_SUCCESS;
}
