/**
 * @file proto_test.h
 * @brief Function prototypes for use with rdtest.c
 * @author Tomas E. Tecce
 */
void load_tables_header(char *);
void load_rscale_tables(double, char *, char *);
double get_data_from_table(double ****, double ****, double, double, 
			   double, double, double);
void free_rscale_tables();

