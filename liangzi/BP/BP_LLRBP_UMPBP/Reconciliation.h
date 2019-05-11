#include <fstream>
using namespace std;

void calNumOfRL(); 
void getZ();
int BP(double e,int iteraM);
int LLRBP(double e,int iteraM);
int UMPBP(double e,int iteraM);
void readH(const char * address);
void randomX();
void randomY(double e);
int Test();
void PostProcess(int * SiftedKey);
double calF(double e);
void concludeDegreeCoefficient(); 
void createXY(int num,double e); 
void readX(int row,const char *address);
void readY(int row,const char *address);
int getN(); 

