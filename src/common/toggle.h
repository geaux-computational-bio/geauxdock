#ifndef  TOGGLE_H
#define  TOGGLE_H


// output toggle
#define IS_H5DUMP 0

// move ligand off the prt center
#define IS_AWAY 0

// control the purterbation be 2.0f
#define IS_CONTROL_MOVE 0

// calculate the total energy using linear combination
#define IS_LINEAR 1

// use bayesian force field
#define IS_BAYE 0

// robust check at run time
#define IS_ASSERT 1

// performance instrument using PAPI library
#define IS_PAPI 0


// remove irregularity in the kernel for performance modeling
#ifndef IS_NO_BRANCH
#define IS_NO_BRANCH 0
#endif

#ifndef CALC_PRT
#define CALC_PRT 1
#endif

#ifndef CALC_KDE
#define CALC_KDE 1
#endif

#ifndef CALC_MCS
#define CALC_MCS 1
#endif

#ifndef CALC_DST
#define CALC_DST 1
#endif





#if IS_NO_BRANCH == 1
#define OROR1 || 1
#else
#define OROR1
#endif





#endif // TOGGLE_H

