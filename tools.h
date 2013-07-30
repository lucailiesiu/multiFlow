#ifndef TOOLS_H
#define TOOLS_H

#define PI ((double)3.1415926535897932385)
#define HALFPI ((double)1.5707963267948966192)
#define TWOPI ((double)6.2831853071795864769)
#define FOURPI ((double)12.566370614359172954)
#define ABS(x) ((x) < 0.0? -(x) : (x))
#define DELTAN ((double)0.001)
#define FIELD_SCALE ((double)1.00)
#define INITIAL_SCALE ((double)1.000)
#define MINSAVE ((int) 5)
#define MAXSAVE ((int) 10000)
#define CLASSICEFOLDS ((int) 60)
#define CE ((double) 0.2703625)
#define C ((double) 0.08145)
#define PARSE ((int) 50)
#define PRINTPARSE ((int) 1000)
#define PARSEFIXED ((int) 500)

enum solutionType {FIXED, ZEROFIXED, NONTRIVIAL, INSUF, SLOWMOVING, NOTFOUND};

#endif
