//radio.h
//misc functions for my radio project

#ifndef RADIO_H
#define RADIO_H


class Radio {

 public:
  Radio(Cosmology *c1, Astrophysics *a1);
  ~Radio();

 protected:

  Cosmology *c;
  Astrophysics *a;
  
};

double dndzdo(double z, Cosmology *c);
double dndo(double z1, double z2, Cosmology *c);
double setDnDo(double z, Cosmology *c1, int iflag);
double getDnDo(double z);

#endif
