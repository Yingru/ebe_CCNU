#ifndef HIJINGINI
#define HIJINGINI

typedef const double CD;
double fws(double r);  //woods-saxon distribution
void wsxy(double &x, double &y);// woods saxon sampling
void generate_xy(double &x, double &y, CD &E, CD &PX, CD &PY, CD &PZ);

#endif
