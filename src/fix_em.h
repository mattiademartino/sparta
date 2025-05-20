/* -*- c++ -*- ----------------------------------------------------------
   SPARTA - Stochastic PArallel Rarefied-gas Time-accurate Analyzer
   http://sparta.sandia.gov
   Steve Plimpton, sjplimp@gmail.com, Michael Gallis, magalli@sandia.gov
   Sandia National Laboratories

   Copyright (2014) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   EM Field fix for SPARTA - A customizable EM field implementation
   Mattia Morganti, m.morganti@em.uni-frankfurt.de
   Institute for Atmospheric Physics, Goethe University Frankfurt
------------------------------------------------------------------------- */
// fix_em.h
#ifdef FIX_CLASS
FixStyle(em,FixEM)
#else
#ifndef SPARTA_FIX_EM_H
#define SPARTA_FIX_EM_H

#include "fix.h"

namespace SPARTA_NS {

class FixEM : public Fix {
 public:
  typedef void (*FieldFunction)(double x, double y, double z, double t, double *fvec);

  enum FieldType { NONE, UNIFORM, CUSTOM, GRID };

  FixEM(class SPARTA *, int, char **);
  ~FixEM();
  int setmask();
  void init();
  void setup();
  void start_of_step();
  void end_of_step();
  void pre_update();
  
  void apply_lorentz_force(int, double *, double *, double, double);
  
  void set_electric_field(double, double, double);
  void set_magnetic_field(double, double, double);
  void set_electric_field_function(FieldFunction);
  void set_magnetic_field_function(FieldFunction);
  void update_custom(int, double, double, double, double *);

 private:
  int field_type;
  double *efield;
  double *bfield;
  FieldFunction efunc;
  FieldFunction bfunc;
  
  int nglocal;
  double **egrid;
  double **bgrid;
  
  double e_charge;    // Elementary charge
  double k_boltz;     // Boltzmann constant
  double mvv2e;       // Mass-velocity conversion factor
  
  double dtfix;       // Time interval for field updates
  
  bool species_with_charge;
  
  double em_time;
  int em_count;
  
  void init_physics();
  void check_species();
  void allocate();
  void compute_em_grid();
  void read_efield(const char *);
  void read_bfield(const char *);
};

}

#endif /* SPARTA_FIX_EM_H */
#endif
