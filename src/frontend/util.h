#ifndef  UTIL_H
#define  UTIL_H


#include <list>
#include <string>


#include <geauxdock.h>

using namespace std;


void Usage (char *);
void Banner ();
void TraceBanner ();
void ParseArguments (int argc, char **argv, McPara *, ExchgPara *, InputFiles *);

void InitLigCoord (Ligand *, const ComplexSize);

void SetPocketCenter (Protein *, const float, const float, const float, const int);
void SetTemperature (Temp *, ExchgPara *);
void SetReplica (ReplicaMC *, const ComplexSize);
void SetMcLog (McLog *);

void DumpRecord (const Record *, const int, const char*);


#endif // UTIL_H

