#include "EBMS.h"
#include "matrix.h"
#include "mpi.h"
#include "runtime_parameters.h"
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<time.h>

int main(int argc, char **argv){
  int nprocs;                          /* total number of processors */
  int mype;                            /* my rank */
  int nt;                              /* num tracking procs */
  int nm;                              /* num memory procs */
  int nb;                              /* num energy bands */
  long npg;
  long lsize,lsizeb;                   /* local size for x-section data (b means "bytes") */
  long gsize,gsizeb;                    /* global size of x-section data in bytes (for consistency) */
  long npl;                            /* local particle count */
  int nranks;                          /* input number of mpi procs -- used for consistency checking */
  float tracking_rate;
  Range *ranges;                       /* max and min energy for each band */
  float **scattering_matrix;           /* intra-group scattering probability */
  Particle *p;                         /* local list of particles */
  float *xdata;
  Boolean tracking_proc = FALSE;       /* am I tracking proc? */
  char sm_file[128];                   /* path to file containing scattering matrix */
  char *param_file;                    /* input parameter file */
  int i;                               /* looping var */
  int size_sm = 0;
  Boolean use_file = FALSE;
  int nonblocking = 0;
  Boolean use_nonblocking =  FALSE;
  MPI_Comm  intracomm, intercomm;
  
  /* initialize all ranks */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
  
  //  srand(time(NULL));
  
  if (mype == MASTER_PE){
      if (argc != 2){
	  printf("must enter input parameter file\n");
	  MPI_Abort(MPI_COMM_WORLD,1);
      }
      param_file = argv[1];
      read_params(param_file);
      list_params();
      assert(get_param("nranks", &nranks) != -1);
      assert(get_param("nt", &nt) != -1);
      assert(get_param("nm", &nm) != -1);
      assert(get_param("nb", &nb) != -1);
      assert(get_param("lsizeb", &lsizeb) != -1);
      assert(get_param("gsizeb", &gsizeb) != -1);
      assert(get_param("npl", &npl) != -1);
      assert(get_param("npg", &npg) != -1);
      assert(get_param("tracking_rate", &tracking_rate) != -1);
      // logical exclusive or operator
      assert((get_param("size_sm", &size_sm) != -1) != (get_param("sm_file", sm_file) != -1));
      printf("nt:%d nm:%d nb:%d lsizeb:%ld npl:%ld scattering_file:%s\n",
	     nt,nm,nb,lsizeb,npl,sm_file);
      /* do some consistency checks */
      if (size_sm == 0) use_file = TRUE;
      get_param("nonblocking", &nonblocking);
      if (nonblocking != 0) use_nonblocking = TRUE;

      else assert(size_sm == nb);
      assert(nranks == nprocs);   /* was mpi started with correct nprocs for problem? */
//      assert(npg/nt == npl);           
//      assert(gsizeb/nm == lsizeb);
      assert (nm + nt == nprocs);
  }
  
  /* send non-redundant input parameters to all procs */
  MPI_Bcast(&tracking_rate, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nt, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&use_nonblocking, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nm, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nb, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&lsizeb, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(&npl, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  
  /* create MPI communicators for disjoint memory and tracking proc sets */
  if (mype < nt) tracking_proc = TRUE;     /* tracking procs are first nt ranks */
  /* create MPI intra- and inter-communicators */
  MPI_Comm_split(MPI_COMM_WORLD,tracking_proc,mype,&intracomm);
  if (tracking_proc)
    MPI_Intercomm_create(intracomm,0,MPI_COMM_WORLD,nt,0,&intercomm);
  else
    MPI_Intercomm_create(intracomm,0,MPI_COMM_WORLD,0,0,&intercomm);

  assert (lsizeb % sizeof(float) == 0); /* test in converting lsizeb to lsize */
  lsize = lsizeb/sizeof(float);         /* do conversion */

  /*allocate space for x-section bands on all procs. On tracking procs
    we need at minimum space for the lasrgest band. Since all bands are
    forced to be equal in this miniapp all allocations are are size lsize */
  assert ( (xdata = (float *) malloc(lsize*sizeof(float))) != NULL);
  
  /* main tracking/data serving loop */
  if (!tracking_proc){
    int ping;
    MPI_Status stat;
    int ncomplete = 0; /* number of finished tracking procs */
    Boolean finished = FALSE;
    MPI_Request *sreqs;
    if (use_nonblocking) sreqs = (MPI_Request *) malloc(nt*sizeof(MPI_Request));

    init_xdata(xdata,lsize);

    while (!finished){
      MPI_Recv(&ping,1,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,
	       intercomm,&stat);
      if (ping == TERMINATION_VALUE){
	++ncomplete;
	printf("memory proc: %d recived -1 from %d\n", mype, stat.MPI_SOURCE);
	if (ncomplete == nt) finished=TRUE;
      }
      else{
        if (use_nonblocking)
          MPI_Isend(xdata,lsize,MPI_FLOAT,stat.MPI_SOURCE,0,intercomm,&sreqs[stat.MPI_SOURCE]);
        else
          MPI_Send(xdata,lsize,MPI_FLOAT,stat.MPI_SOURCE,0,
                   intercomm);
      }
    }
    printf("memory proc %d shutting down ...\n",mype);
  }
  else{              /* if tracking proc */
    int i,j,k;
    double ran_val;
    double ts,tf, ts_comm, tf_comm, comm_time;
    comm_time = 0.0;
    long *n_alive = (long *) malloc(nm*sizeof(long));
    MPI_Status stat;
    
    int ping = 0;                  /* used to communicate intention to server procs */
    int n_passes = 0;              /* number of scattering matrix sweeps */
    long total_alive;              /* total particles not absorbed */
    float probability;             /* used to assign outscatter group */

    struct timespec pause, rem;    /* sleep time to imitate particles */

    /* intialize scatterin matrix. band 0 is highest energy, nb-1 lowest */
    scattering_matrix = matrix(0,nb-1,0,nb-1);

    /* read scattering matrix on master and broadcast to tracking procs only */
    if (mype == MASTER_PE){
      if (use_file == TRUE){
	printf("reading scattering matrix from file....\n");
	read_scattering_matrix(sm_file,scattering_matrix,nb);
      }
      else {
	printf("creating scattering matrix....\n");
	create_scattering_matrix(scattering_matrix,nb);     
      }
    }

    MPI_Bcast(&scattering_matrix[0][0],nb*nb,MPI_FLOAT,MASTER_PE,intracomm);

    /* set energy ranges associated with each band (min..max). This
       isn't really needed in this communication kernel but is included
       for conceptual purposes in case we choose to extend kernel */
    assert( (ranges =(Range *) malloc(nm*sizeof(Range))) != NULL);
    set_energy_ranges(ranges,nm);

    /* allocate local list of particles and initialize all to highest
       energy band (ie band 0). Exact energies are not important for this
       kernel, only band.
    */
    assert( (p = (Particle *) malloc(npl*sizeof(Particle))) != NULL);
    init_particles(p,npl,mype);

    /* need to count how many alive in each band so
       that we don't load memory on subsequent passes
       when there are none alive. On first pass all neutrons
       are in band 0.
    */
    n_alive[0] = npl;
    for (i=1;i<nb;++i) n_alive[i] = 0.;

    /* loop over bands starting with highest energy */
    total_alive = npl;

    ts = MPI_Wtime();
    pause.tv_nsec = (long) (tracking_rate*1.e6);    /* tracking rate is in milliseconds */	 
    pause.tv_sec = (time_t) 0;
    int interaction_num = 0;

    /*temp check variables */
    double start_time,end_time;
    int particle_deaths=0;

    while (total_alive > 0){
      ++n_passes;
      for (k=0;k<nb;++k){
	int server_proc = k;
	ts_comm = MPI_Wtime();
	MPI_Send(&ping,1,MPI_INTEGER,server_proc, 0,
		 intercomm);
	
	MPI_Recv(xdata,lsize,MPI_FLOAT,
		 server_proc, MPI_ANY_TAG,
		 intercomm, &stat);
	tf_comm = MPI_Wtime();
	comm_time += (tf_comm - ts_comm);
	/* Temporary check for hypothesis, only band 0 data */
	if (k==1) printf("deaths %d effective absorption rate: %lf time to request: %lf\n",particle_deaths,((double) particle_deaths)/((double) npl),end_time-start_time);

	if (k==0) start_time = MPI_Wtime();
	for (i=0;i<npl;++i){
	    while ( p[i].band == k && !p[i].absorbed){
		ran_val = rn();
		interaction_num++;
		if (ran_val <= ABSORPTION_THRESHOLD){
		    assert (nanosleep(&pause,&rem) == 0); 
		    p[i].absorbed = TRUE;
		    // temporary check
		    if (k ==0) particle_deaths++;
		    --n_alive[k];
		}
		else{
		    ran_val = rn();
		    probability = scattering_matrix[k][0];
		    j = 0;
		    while (ran_val >= probability){
			++j;
			probability = scattering_matrix[k][j];
		    }
		    p[i].band = j;
		    --n_alive[k];++n_alive[j];
		    
		}
	    }
	}
	
	total_alive = tot(n_alive,nm);
	if(k==0) end_time = MPI_Wtime();
      }
    }

    tf=MPI_Wtime();
    printf("simulation terminated: proc:%d  passes:%d time(s):%f comm_time(s):%f, \n", 
	   mype, n_passes,tf-ts,comm_time);
    
    for (i=0;i<nm;++i){
      ping=-1; /*shutdown signal */
      MPI_Send(&ping,1,MPI_INT,i,99,intercomm);
    }
    FILE* file;
    double max_comm_time,min_tracking_time,max_tracking_time, min_comm_time;
    double total_time = tf -ts;
    MPI_Reduce(&comm_time,&max_comm_time,1,MPI_DOUBLE,MPI_MAX,MASTER_PE,intracomm);
    MPI_Reduce(&comm_time,&min_comm_time,1,MPI_DOUBLE,MPI_MIN,MASTER_PE,intracomm);
    MPI_Reduce(&total_time,&min_tracking_time,1,MPI_DOUBLE,MPI_MAX,MASTER_PE,intracomm);
    MPI_Reduce(&total_time,&max_tracking_time,1,MPI_DOUBLE,MPI_MIN,MASTER_PE,intracomm);
    if (mype == MASTER_PE){
      file = fopen("summary.out","a");
      fprintf(file,"%d \t%d \t%ld \t%ld \t%lf \t%lf \t%lf \t%lf \t%lf \t%d \n",nt,nm,npl,lsizeb,max_comm_time,min_comm_time,npl/min_tracking_time,npl/max_tracking_time,ABSORPTION_THRESHOLD,INT_PER_NEUTRON);
      fclose(file);
    }
  }
  MPI_Finalize();
  
}

/*-----------------------------------------*/
static void read_scattering_matrix(char *sm_file, float **scattering_matrix, int nb){
  int i,j, n, m;
  float val;
  FILE * file;
  int newline;
  file=fopen(sm_file,"r");
  fscanf(file,"%d",&n);
  fscanf(file,"%d",&m);
  assert (n == nb);
  assert (n == m);
  
  for(i=0;i<n;i++) {
    for(j=0;j<m;j++) {
      fscanf(file,"%f",&val);
      scattering_matrix[i][j] = val;
      printf("%f ", scattering_matrix[i][j]);
    }
    printf("%\n");
  }
    
  fclose(file);
}
/* END MAIN */

 
/*---------------------------------------------------
  - void init_xdata(inout xdata, in nl)
  intialiaze cross section data. for purposes
  of this kernel app the values are arbitrary, only
  size of data structure matters. 
  ------------------------------------------------------*/
static void init_xdata(float *xdata, int nl){
  int i;
    for (i = 0; i < nl; ++i)
      xdata[i] = 1.0;
}


static void init_particles(Particle *p, long npl, int mype){
  long i;
  for (i=0;i<npl;++i){
    p[i].energy = 2.0;
    p[i].absorbed = FALSE;
    p[i].band = 0;
    p[i].proc = mype;
  }
}

static void set_energy_ranges(Range *energy_ranges, int nm){
  int i;
  float delta = 2.0/nm;
  energy_ranges[0].lower=0.0;
  energy_ranges[0].upper=delta;
  for (i=1;i<nm;++i){
    energy_ranges[i].lower = energy_ranges[i-1].upper;
    energy_ranges[i].upper = energy_ranges[i].lower + delta;
  }   
}

static long tot(long f[], long n){
  long sum = 0;
  int i;
  for (i=0;i<n;++i)
    sum += f[i];
  return sum;
}

static void create_scattering_matrix(float **scattering_matrix, int nb){
  int i,j;

  //create pdf                                                                           
  for (i = 0; i < nb; ++i)
    for (j = i; j < nb; ++j)
      scattering_matrix[i][j] = 1./(nb-i);

  //convert to cdf for mini-app                                                          
  for (i = 0; i < nb; ++i){
    for (j = 1; j < nb; ++j){
      scattering_matrix[i][j] += scattering_matrix[i][j-1];
      printf("%f ", scattering_matrix[i][j]);
    }
    printf("\n");
    scattering_matrix[i][nb-1] = 1.0;
  }
  for (i = 0; i < nb; ++i){
    for (j = 0; j < nb; ++j){
      printf("%f ", scattering_matrix[i][j]);
    }
    printf("\n");
  }
}
