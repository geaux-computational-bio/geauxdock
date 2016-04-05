#ifndef  UTIL_PRINT_H
#define  UTIL_PRINT_H


#include <cstring>
#include <string>

#include <size.h>
#include <geauxdock.h>


void PrintTrack (Record *,  int,  int,  int,  int);


void PrintMoveVector (const float*, const int);
void PrintMoveRecord (const Record *, const int, const int, const int, const int, const int);
void PrintCsv (const float *, const int, const int, const int);

void PrintRecord (Record *,  int,  int,  int,  int,  int);
void PrintRepRecord (const Record *, const int, const int, const int, const int, const int, const int);
void PrintRepRecord2 (Record *,  ComplexSize,  int,  int,  int,  int,  int,  int);
void PrintLigCoord (const Ligand *, const int);
void PrintLigand (const Ligand *);
void PrintProtein (const Protein *);
void PrintDataSize (const Complex *);
void PrintSummary (const Complex *);


#endif // UTIL_PRINT_H


