#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>
#include <list>
#include <sys/stat.h>


#include <size.h>
#include <toggle.h>
#include <geauxdock.h>
#include <util.h>
#include <load.h>

#include <yeah/c/mkdir.h>




using namespace std;


void
Usage (char *bin)
{
  fprintf (stderr, "usage: %s [options]\n", bin);
  exit (1);
}






void
Banner ()
{
  cout << "------------------------------------------------------------" << endl
       << "                         GeauxDock                          " << endl
       << "                        version 0.1                         " << endl << endl
       << "   GPU-accelerated mixed-resolution ligand docking using    " << endl
       << "                ReplicaMC Exchange Monte Carlo                " << endl
       << "------------------------------------------------------------" << endl << endl;
}




void
ParseArguments (int argc, char **argv, McPara * mcpara, ExchgPara * exchgpara,
		InputFiles * inputfiles)
{
  float ts, rs;

  // essential default parameters
#if IS_H5DUMP == 1
  char mydir[MAXSTRINGLENG] = "output_default";
#endif
  mcpara->steps_total = STEPS_TOTAL;
  mcpara->steps_per_dump = STEPS_PER_DUMP;
#if IS_BAYE == 1
  inputfiles->norpara_file.path_a = "baye_nor_a";
  inputfiles->norpara_file.path_b = "baye_nor_b";
#elif IS_BAYE == 0
  inputfiles->norpara_file.path_a = "08_nor_a";
  inputfiles->norpara_file.path_b = "08_nor_b";
#endif


#if 1
  // default parameters for debug convenience

  std::string basedir ("../data/");
  inputfiles->lig_list = basedir + "astex/ligands/ligs.txt";
  inputfiles->prt_list = basedir + "astex/proteins/prt.txt";
  inputfiles->weight_file.path = basedir + "parameters/08ff_opt";
  inputfiles->enepara_file.path = basedir + "parameters/paras";
  inputfiles->norpara_file.path_a = basedir + "parameters/08_nor_a";
  inputfiles->norpara_file.path_b = basedir + "parameters/08_nor_b";
  exchgpara->floor_temp = 0.044f;
  exchgpara->ceiling_temp = 0.036f;
  exchgpara->num_temp = 1;
  ts = 0.02f;
  rs = 0.08f;
#endif  




  for (int i = 0; i < argc; i++) {

    // lig list
    if (!strcmp (argv[i], "-ll") && i < argc) {
      inputfiles->lig_list = argv[i + 1];
    }
    // prt list
    if (!strcmp (argv[i], "-lp") && i < argc) {
      inputfiles->prt_list = argv[i + 1];
    }
    // parameter files
    if (!strcmp (argv[i], "-opt") && i < argc) {
      inputfiles->weight_file.path = argv[i + 1];
    }
    if (!strcmp (argv[i], "-na") && i < argc) {
      inputfiles->norpara_file.path_a = argv[i + 1];
    }
    if (!strcmp (argv[i], "-nb") && i < argc) {
      inputfiles->norpara_file.path_b = argv[i + 1];
    }
    if (!strcmp (argv[i], "-para") && i < argc) {
      inputfiles->enepara_file.path = argv[i + 1];
    }


    // MC steps
    if (!strcmp (argv[i], "-ns") && i < argc) {
      mcpara->steps_total = atoi (argv[i + 1]);
    }


    // temperatures
    if (!strcmp (argv[i], "-floor_temp") && i < argc) {
      exchgpara->floor_temp = atof (argv[i + 1]);
    }
    if (!strcmp (argv[i], "-ceiling_temp") && i < argc) {
      exchgpara->ceiling_temp = atof (argv[i + 1]);
    }
    if (!strcmp (argv[i], "-nt") && i < argc) {
      int num_temp = atoi (argv[i + 1]);
      if (num_temp == 0) {
	cout << "number of temperature cannot set to be zero" << endl;
	cout << "docking exiting ..." << endl;
	exit (1);
      }
      else if (num_temp == 1) {
	exchgpara->num_temp = num_temp;
      }
      else {
	if ((num_temp <= MAX_TMP) && (num_temp > 1))
	  exchgpara->num_temp = num_temp;
	else {
	  cout << "setting number of temperatures exceeds MAX_TMP" << endl;
	  cout << "try modifying MAX_TMP in size.h and compile again" << endl;
	  cout << "docking exiting ..." << endl;
	  exit (1);
	}
      }
    }
    
    // trace file
    if (!strcmp (argv[i], "-tr") && i < argc) {
      inputfiles->trace_file.path = argv[i + 1];
    }

    // move scale
    if (!strcmp (argv[i], "-t") && i < argc) {
      ts = atof (argv[i + 1]);
    }
    if (!strcmp (argv[i], "-r") && i < argc) {
      rs = atof (argv[i + 1]);
    }
  }

  /*
  if (!protein_opt) {
    cout << "Provide target protein structure" << endl;
    exit (EXIT_FAILURE);
  }

  if (!compounds_opt) {
    cout << "Provide compound library in SD format" << endl;
    exit (EXIT_FAILURE);
  }

  if (!lhm_opt) {
    cout << "Provide LHM potentials" << endl;
    exit (EXIT_FAILURE);
  }
  */

  mcpara->move_scale[0] = ts;
  mcpara->move_scale[1] = ts;
  mcpara->move_scale[2] = ts;
  mcpara->move_scale[3] = rs;
  mcpara->move_scale[4] = rs;
  mcpara->move_scale[5] = rs;


#define MAX_LINE 1024
  ifstream fnl, fnl2;
  int conf_lig = 0;
  char line [MAX_LINE];

  // get lig numbers from the list
  fnl.open (inputfiles->lig_list.c_str());
  if (fnl.fail ()) {
    std::cout << "fail to open ligand list" << std::endl;
    exit (1);
  }
  while (fnl) {
    fnl.getline (line, MAX_LINE);
    char delim[] = " ";
    char *token0;
    token0 = strtok (line, delim);
    if (token0 != NULL) {
      conf_lig++;
    }
  }
  fnl.close();


  // allocate ligand strings
  inputfiles->lig_files_num = conf_lig;
  inputfiles->lig_files = new LigandFile[conf_lig];
  inputfiles->lhm_files = new LhmFile[conf_lig];



  // loading ligand files from the list
  fnl2.open (inputfiles->lig_list.c_str());
  if (fnl2.fail ()) {
    std::cout << "fail to open ligand list" << std::endl;
    exit (1);
  }
  int counter = 0;
  while (fnl2) {
    fnl2.getline (line, MAX_LINE);
    char delim[] = " ";
    char *token0;
    char *token1;
    char *token2;
    token0 = strtok (line, delim);
    token1 = strtok (NULL, delim);
    token2 = strtok (NULL, delim);
    if (token0 != NULL) {
      std::string s0(token0);
      std::string s1(token1);
      std::string s2(token2);
      inputfiles->lhm_files[counter].ligand_id = s0;
      inputfiles->lig_files[counter].path = s1;
      inputfiles->lhm_files[counter].path = s2;
      inputfiles->lig_files[counter].molid = "MOLID";
      //std::cout << s0 << "end" << std::endl;
      counter++;
    }
  }
  fnl2.close();



  ifstream fnp;
  fnp.open (inputfiles->prt_list.c_str());
  if (fnp.fail ()) {
    std::cout << "fail to open protein list" << std::endl;
    exit (1);
  }
  while (fnp) {
    fnp.getline (line, MAX_LINE);
    char delim[] = " ";
    char *token0;
    char *token1;
    token0 = strtok (line, delim);
    token1 = strtok (NULL, delim);
    if (token0 != NULL) {
      std::string s0(token0);
      std::string s1(token1);
      inputfiles->prt_file.id = s0;
      inputfiles->prt_file.path = s1;
    }
  }
  fnp.close();









#if IS_H5DUMP == 1
  // create output directory

  // obtain the time tag
  char mystime[MAXSTRINGLENG];
  time_t mytime = time (NULL);
  struct tm *mylocaltime;
  mylocaltime = localtime (&mytime);
  strftime (mystime, MAXSTRINGLENG, "%Y%m%d_%H%M%S", mylocaltime);

  // name the output directory using time tag
  if (strcmp (mydir, "output_default") == 0) {
    sprintf (mydir, "output_%s", mystime);
  }

  strcpy (mcpara->outputdir, mydir);
  MakeDir (mydir);

  // prefix of the file name
  const char h5file[MAXSTRINGLENG] = "a";
  strcpy (mcpara->outputfile, h5file);
#endif


}



void
SetPocketCenter (Protein * prt, const float cx, const float cy, const float cz, const int n_prt)
{
 for (int i = 0; i < n_prt; ++i) {
    Protein *dst = &prt[i];
    dst->pocket_center[0] = cx;
    dst->pocket_center[1] = cy;
    dst->pocket_center[2] = cz;
 }
}



void
SetTemperature (Temp * temp, ExchgPara * exchgpara)
{
  int num_temp = exchgpara->num_temp;
  float floor_temp = exchgpara->floor_temp;
  float ceiling_temp = exchgpara->ceiling_temp;

  if (num_temp == 1) {
    float my_temp = floor_temp;
    float beta = 1.0f / (BOLTZMANN_CONST * my_temp);

    for (int i = 0; i < num_temp; i++) {
      //printf ("temp # %d\t\t\t%f\n", i, my_temp);
      temp[i].order = i;
      temp[i].minus_beta = 0.0f - beta;
    }
  }
  else {
    const float temp_ratio = exp (log (ceiling_temp / floor_temp) / (float) (num_temp - 1));
    float a = floor_temp;
    for (int i = 0; i < num_temp; i++) {
      float my_temp = a;
      float my_beta = 1.0f / (BOLTZMANN_CONST * my_temp);
      //printf ("temp # %d\t\t\t%f\n", i, my_temp);

      temp[i].order = i;
      temp[i].minus_beta = 0.0f - my_beta;

      a *= temp_ratio;
    }
  }


  // for (int i = 0; i < num_temp; i++) {
  //   temp[i].t = floor_temp;
  //   temp[i].minus_beta = -1.0f / temp[i].t;
  //   temp[i].order = i;
  // }
}




// replica[n_rep]
// replica[n_prt][n_tmp][n_lig]
void
SetReplica (ReplicaMC * replica, const ComplexSize complexsize)
{
  const int n_lig = complexsize.n_lig;
  const int n_prt = complexsize.n_prt;
  const int n_tmp = complexsize.n_tmp;

  for (int p = 0; p < n_prt; ++p) {
    for (int t = 0; t < n_tmp; ++t) {
      for (int l = 0; l < n_lig; ++l) {
	const int flatten_addr = n_tmp * n_lig * p + n_lig * t + l;
	replica[flatten_addr].idx_rep = flatten_addr;


	replica[flatten_addr].idx_prt = p;
	replica[flatten_addr].idx_tmp = t;
	replica[flatten_addr].idx_lig = l;

        /*
	// set identical replicas for test purpose
	replica[flatten_addr].idx_prt = 0;
	replica[flatten_addr].idx_tmp = 0;
	replica[flatten_addr].idx_lig = 0;
        */


	for (int aa = 0; aa < 6; ++aa)
	  replica[flatten_addr].movematrix[aa] = 0.0f;
	for (int aa = 0; aa < MAXWEI; ++aa)
	  replica[flatten_addr].energy[aa] = 0.0f;

	replica[flatten_addr].is_accept = 0;
      }
    }
  }

}


void
SetMcLog (McLog * mclog)
{
  mclog->ac_mc = 0;  // MC acceptance counter
  for (int i = 0; i < MAX_REP; ++i)
    mclog->acs_mc[i] = 0;   // MC acceptance counter for all replicas
}



/*
float MyRand()
{
	//float randdd = (float) rand () / RAND_MAX;
	float randdd = 0.002f;
	return randdd;
}
*/

