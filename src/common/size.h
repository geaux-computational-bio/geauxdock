
/*
==============================================================================================
     __________________ ____ ___              .___             __      _________   _____   
    /  _____/\______   \    |   \           __| _/____   ____ |  | __ /   _____/  /     \  
   /   \  ___ |     ___/    |   /  ______  / __ |/  _ \_/ ___\|  |/ / \_____  \  /  \ /  \ 
   \    \_\  \|    |   |    |  /  /_____/ / /_/ (  <_> )  \___|    <  /        \/    Y    \
    \______  /|____|   |______/           \____ |\____/ \___  >__|_ \/_______  /\____|__  /
           \/                                  \/           \/     \/        \/         \/ 

      GPU-accelerated hybrid-resolution ligand docking using ReplicaMC Exchange Monte Carlo

==============================================================================================
*/


#ifndef __SIZE_H_
#define __SIZE_H_


#define MAX_STR_LENG 256



#define NAMEMAX 64
/* max lig aton name */

#define MAXPRO 2048
/* protein residues */

#define MAXLIG 64
/* ligand heavy atoms */

#define MAX_CONF_PRT  20
/* protein confs */

#define MAX_CONF_LIG 64
/* ligand confs */

#define MAXLIG_NUM 200
/* MAX ligand in one ligand .sdf file */

#define MAX_TMP 1024
//#define MAX_TMP 32
//#define MAX_TMP 1
/* number of temperature replicas */

#define MAX_REP 16384
//#define MAX_REP 4096
//#define MAX_REP 1024
/* max replicas */

#define  MAXSWP 2000
/* max swapping pairs */

#define  MAXLIB 100
/* library cmps */

#define  MAXSDF  500
/* sdf length */

#define MAXTP1 31 // using 32 does not help imporoving performance
/* point types */

#define MAXTP2 24
/* atom types */

#define MAXTP3 50
/* point types (for ele) */

#define MAXTP4 20
/* residue types */

#define MAXFP1 1024
/* smiles */

#define  MAXFP2 168
/* maccs */

// MAXWEI must be less than 16 (can be extended to 32), see CUDA reduction code "e_s"
#define  MAXWEI 10
/* energy terms */

#define MAXKDE 10240
/* kde points */

#define MAX_MCS_ROW 512
/* position restraints */

#define MAX_MCS_COL MAXLIG
/* mcs fields, number of field in a mcs */


#define INITTEMP 10000.0f
/* temperature in the first replica */


#define PI 3.1415926535f

/* value used for return statement in data.C for missing atoms */
#define BADKDE 50

// maximum string length for file names
#define MAXSTRINGLENG 128

// if mcs equal to 123456789.0f, it is invalid
#define MCS_INVALID_COORD 0.01234f

// if mcs equal to 728492, it is invalid
#define CMCC_INVALID_VAL 10000.0f

// boltzman constant
#define BOLTZMANN_CONST 1.0f

// monte carlo steps
#define STEPS_PER_DUMP 200

// steps that no not have to be defined as macro
#define STEPS_TOTAL 200


// signal the MPI client to terminate
#define FINISH_SIGNAL 0XFFFFFF

#endif

