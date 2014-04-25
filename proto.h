/**
 * @file proto.h
 * @brief Function prototypes to include with rdtables.c
 * @author Tomas E. Tecce
 */
void readparameterfile(char *);
void checkforerror(char *, int, char *);

void output_times(void);
void allocate_arrays(void);
void free_arrays(void);

void determine_axes_values(void);


int nrsolver(int, double, double, double, double, double, double, double *);

double compute_deltavir(double);
/*double m_of_rinit(double, double, double, double, double, double,
  double, double, double);*/
double m_of_rinit(double, double, double, double, double, double,
		  double, double);  /* Version for rditer2 */
double rdaprox(double, double, double, double, double);
double mass_nfw(double, double, double);

void test_mode(void);

void load_concentration_tables(char *);
void get_halo_concentration(double, double, double *);
void find_table_limits(void);
void concentration_range(double, double *, double *);
void free_concentration_tables(void);
void get_mass_from_concentration(double, double, double *);

int generate_rscale_table(int, double, char *, char *);

void rditer(double, double, double, double, double, double, double *, 
	    double *, double *);

/*int dzbrac(double *, double *, double, double, double, 
  double, double, double, double, double);*/
int dzbrac(double *, double *, double, double, double, 
	   double, double, double, double); /* Version for rditer2 */
