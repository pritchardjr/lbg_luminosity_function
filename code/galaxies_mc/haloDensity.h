double urhoNFW(double r, double c, double rvir);
double uNFW(double k, double c, double rs);
double concentration(double m, double z);
double findRs(double m, double conc, double z, Cosmology *c);

double PS1hdd(double k, double z, Cosmology *c);
void setPS1hddIntegrand(double m, double y[], double deriv[], double k1, 
			double z1, Cosmology *c1, int flag);
void ps1hddIntegrand(double m, double y[], double deriv[]);
double PS2hdd(double integral, double renorm, double Plin);
double PS2hddFull(double k, double z, double renorm, Cosmology *c);
double PS2hddInt(double k, double z, Cosmology *c);
void setPS2hddIntegrand(double m, double y[], double deriv[], double k1, 
			double z1, Cosmology *c1, int flag);
void ps2hddIntegrand(double m, double y[], double deriv[]);
double RenormInt(double z, Cosmology *c);
void setRenormIntegrand(double m, double y[], double deriv[], double z1,
			Cosmology *c1, int flag);
void renormIntegrand(double m, double y[], double deriv[]);

double xidd(double r, double z, Cosmology *c, int flag);
void setCorrInt(double kp, double y[], double deriv[], double r1, 
		double z1, Cosmology *c1, int flag);
void corrInt(double kp, double y[], double deriv[]);
double PSddTable(double kp, double z, Cosmology *c, int flag);

double Plinear(double k, double z, Cosmology *c);

///
double haloBias(double mass, double z, Cosmology *c);
double meanHaloBias(double z, Cosmology *c);
void setHaloBiasIntegrand(double m, double y[], double deriv[], 
			  double z1, Cosmology *c1, int flag);
void haloBiasIntegrand(double m, double y[], double deriv[]);

double meanHaloBiasNW(double z, Cosmology *c, double mass= -1.0);
void setHaloBiasIntegrandNW(double m, double y[], double deriv[], 
			  double z1, Cosmology *c1, int flag);
void haloBiasIntegrandNW(double m, double y[], double deriv[]);
double biasTinker(double mass, double z, Cosmology *cosm);

double meanHaloNumber(double z, Cosmology *c, double mass= -1.0);
void setHaloNumberIntegrand(double m, double y[], double deriv[], 
			  double z1, Cosmology *c1, int flag);
void haloNumberIntegrand(double m, double y[], double deriv[]);
