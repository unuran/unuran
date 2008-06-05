
struct siv{
double *ui;//[g+1];
double *zi;//[g+1];
double xi;
double cdfi;
};

struct genobject{
 struct siv *iv;//[maxint+1] for setup; for sampling [ni+1]
 int g; //degree of polynomial
 int ni;//number of sub intervals
 int *gt;//[C]   guidetable
 int C; //size of guide table
 double umax;
  
};

struct genobject *pinvsetup( double (*f)(double x),int g, double uerror, double x0, int asearch, int bsearch,double a, double b);
                             

struct genobject *setup(double (*f)(double x),int g, double a, double b, double hh, double uerror);




double lobato5(double x, double h, double fx, double *fxph, double (*f)(double x));

double nint_12(double a,double b,double *res_relerror, double (*f)(double x));

double nint_monoton_dens(double a,double b,double step,double crit, double (*f)(double x));

double evalnewtoninterpol(double u,int g,double ui[],double zi[]);

int newtoninterpol(double x0, double h,int g,double ui[],double zi[],double *x,double (*f)(double x));

int tstpt(int g,double ui[],double utest[]);

double maxerrornewton(int g,double ui[],double zi[],double x0,double xval[],double (*f)(double x));

struct genobject *init_genobject(int g,int maxint);

int free_genobject(struct genobject *p);

double searchborder(double x0, double step,double border,double (*f)(double x));

double tail(double x, double d,double (*f)(double x));

double cut(double w,double dw, double crit,double (*f)(double x));

double quantile(double u, struct genobject *geno);




