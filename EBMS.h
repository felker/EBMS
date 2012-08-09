#define ABSORPTION_THRESHOLD .1
#define TERMINATION_VALUE -1
#define rn() (rand())/((double) RAND_MAX)

#define MASTER_PE 0
#define MAX_FILENAME_CHARS = 128
typedef enum boolean{
  FALSE, 
  TRUE
} Boolean;

typedef struct range{
  float lower;
  float upper;
} Range;

typedef struct proc{
  int       lrank;
  int       grank;
  int       type;
} Proc;


typedef struct Particle_type{
  int band;
  double energy;
  Boolean absorbed;
  int proc;
} Particle;

float **matrix(long nrl, long nrh, long ncl, long nch);
static void create_scattering_matrix(float **scattering_matrix, int nb);
static void read_scattering_matrix(char* sm_file, float **scattering_matrix, int nb);
static void init_xdata(float *xdata, int nl);
static void init_particles(Particle *p, long npl, int mype);
static void set_energy_ranges(Range*, int);
static void process_input(int argc, char **argv, int *nt, int *nm, int *nb);
static long tot(long f[], long n);


