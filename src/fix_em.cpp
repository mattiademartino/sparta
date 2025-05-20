// fix_em.cpp
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <cstdio>  // For snprintf
#include "fix_em.h"
#include "update.h"
#include "particle.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "sparta_masks.h"
#include "grid.h"
#include "math_const.h"

#ifndef PRE_UPDATE
#define PRE_UPDATE 2
#endif

using namespace SPARTA_NS;
using namespace MathConst;

#ifdef FIX_CLASS
FixStyle(em,FixEM)
#endif

// helper for formatted error
static void error_all(Error *error, const char *fmt, const char *arg) {
  char msg[128];
  snprintf(msg, 128, fmt, arg);
  error->all(FLERR, msg);
}

FixEM::FixEM(SPARTA *sparta, int narg, char **arg)
  : Fix(sparta, narg, arg),
    field_type(NONE), efield(nullptr), bfield(nullptr),
    efunc(nullptr), bfunc(nullptr),
    nglocal(0), egrid(nullptr), bgrid(nullptr),
    e_charge(0.0), k_boltz(0.0), mvv2e(0.0),
    dtfix(0.0), species_with_charge(false),
    em_time(0.0), em_count(0)
{
  if (narg < 2) error->all(FLERR, "FixEM: too few args");
  init_physics();

  efield = new double[3];
  bfield = new double[3];
  for (int i = 0; i < 3; i++) {
    efield[i] = 0.0;
    bfield[i] = 0.0;
  }

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "efield") == 0) {
      if (iarg + 3 >= narg) error->all(FLERR, "FixEM: not enough efield values");
      efield[0] = std::atof(arg[iarg + 1]);
      efield[1] = std::atof(arg[iarg + 2]);
      efield[2] = std::atof(arg[iarg + 3]);
      if (field_type < UNIFORM) field_type = UNIFORM;
      iarg += 4;
      if (comm->me == 0) {
        char msg[128];
        snprintf(msg, 128, "FixEM: parsed uniform efield = (%g, %g, %g)",
                 efield[0], efield[1], efield[2]);
        error->message(FLERR, msg);
      }
    } else if (strcmp(arg[iarg], "bfield") == 0) {
      if (iarg + 3 >= narg) error->all(FLERR, "FixEM: not enough bfield values");
      bfield[0] = std::atof(arg[iarg + 1]);
      bfield[1] = std::atof(arg[iarg + 2]);
      bfield[2] = std::atof(arg[iarg + 3]);
      if (field_type < UNIFORM) field_type = UNIFORM;
      iarg += 4;
      if (comm->me == 0) {
        char msg[128];
        snprintf(msg, 128, "FixEM: parsed uniform bfield = (%g, %g, %g)",
                 bfield[0], bfield[1], bfield[2]);
        error->message(FLERR, msg);
      }
    } else if (strcmp(arg[iarg], "efile") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "FixEM: missing efile");
      read_efield(arg[iarg + 1]);
      if (field_type < GRID) field_type = GRID;
      iarg += 2;
      if (comm->me == 0) {
        char msg[128];
        snprintf(msg, 128, "FixEM: reading efield from file '%s'", arg[iarg - 1]);
        error->message(FLERR, msg);
      }
    } else if (strcmp(arg[iarg], "bfile") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "FixEM: missing bfile");
      read_bfield(arg[iarg + 1]);
      if (field_type < GRID) field_type = GRID;
      iarg += 2;
      if (comm->me == 0) {
        char msg[128];
        snprintf(msg, 128, "FixEM: reading bfield from file '%s'", arg[iarg - 1]);
        error->message(FLERR, msg);
      }
    } else if (strcmp(arg[iarg], "update") == 0) {
      if (iarg + 1 >= narg) error->all(FLERR, "FixEM: missing update interval");
      dtfix = std::atof(arg[iarg + 1]);
      if (dtfix <= 0.0) error->all(FLERR, "FixEM: update must be positive");
      iarg += 2;
      if (comm->me == 0) {
        char msg[128];
        snprintf(msg, 128, "FixEM: update interval dtfix = %g", dtfix);
        error->message(FLERR, msg);
      }
    } else {
      error_all(error, "FixEM: unknown keyword '%s'", arg[iarg]);
    }
  }
  if (field_type == NONE)
    error->warning(FLERR, "FixEM: no fields specified");
}

FixEM::~FixEM()
{
  delete[] efield;
  delete[] bfield;
  if (egrid) memory->destroy(egrid);
  if (bgrid) memory->destroy(bgrid);
  if (em_count > 0 && comm->me == 0) {
    char msg[128];
    snprintf(msg, 128, "EM: %d calls in %g s (%g s/call)",
             em_count, em_time, em_time / em_count);
    error->message(FLERR, msg);
  }
}

int FixEM::setmask() {
  return PRE_UPDATE;
}

void FixEM::init()
{
  if (comm->me == 0) error->message(FLERR, "FixEM: init() called");
  check_species();
  if (comm->me == 0) {
    for (int i = 0; i < particle->nspecies; i++) {
      char msg[256];
      snprintf(msg, 256, "FixEM: specie %d mass = %g, charge = %g",
               i, particle->species[i].mass, particle->species[i].charge);
      error->message(FLERR, msg);
    }
  }
  if (field_type == GRID) {
    if (!grid) error->all(FLERR, "FixEM: grid undefined");
    allocate();
    compute_em_grid();
  }
}

void FixEM::setup()
{
  if (comm->me == 0) error->message(FLERR, "FixEM: setup() called");
}

void FixEM::start_of_step()
{
  if (comm->me == 0) error->message(FLERR, "FixEM: start_of_step() called");
}

void FixEM::end_of_step()
{
  if (comm->me == 0) error->message(FLERR, "FixEM: end_of_step() called");
}

void FixEM::pre_update()
{
  if (comm->me == 0) error->message(FLERR, "FixEM: pre_update() called");

  if (!species_with_charge) {
    if (comm->me == 0) error->warning(FLERR, "FixEM: no charged species, skipping pre_update");
    return;
  }

  double dt = update->dt;
  if (comm->me == 0) {
    char msg[128];
    snprintf(msg, 128, "FixEM: pre_update with dt = %g and %d local particles",
             dt, particle->nlocal);
    error->message(FLERR, msg);
  }

  for (int i = 0; i < particle->nlocal; i++) {
    Particle::OnePart &p = particle->particles[i];
    // Check if species index is valid
    if (p.ispecies < 0 || p.ispecies >= particle->nspecies) {
      if (comm->me == 0) {
        char msg[128];
        snprintf(msg, 128, "FixEM: invalid species index %d for particle %d", 
                 p.ispecies, i);
        error->warning(FLERR, msg);
      }
      continue;
    }
    double q = particle->species[p.ispecies].charge;
    double m = particle->species[p.ispecies].mass;
    if (std::fabs(q) < 1e-12 || m <= 0.0) {
      if (comm->me == 0) {
        char msg[128];
        snprintf(msg, 128, "FixEM: skip particle %d (q=%g, m=%g)", i, q, m);
        error->message(FLERR, msg);
      }
      continue;
    }

    if (comm->me == 0 && i == 0) {
      char msg[256];
      snprintf(msg, 256,
               "FixEM: applying lorentz force on particle %d: pos=(%g,%g,%g), vel=(%g,%g,%g), q=%g, m=%g",
               i, p.x[0], p.x[1], p.x[2], p.v[0], p.v[1], p.v[2], q, m);
      error->message(FLERR, msg);
    }

    apply_lorentz_force(i, efield, bfield, q, m);

    if (comm->me == 0 && i == 0) {
      Particle::OnePart &p_after = particle->particles[i];
      char msg[256];
      snprintf(msg, 256,
               "FixEM: after lorentz force: vel=(%g,%g,%g)",
               p_after.v[0], p_after.v[1], p_after.v[2]);
      error->message(FLERR, msg);
    }
  }
}

void FixEM::apply_lorentz_force(int ip, double *evec, double *bvec, double q, double m)
{
  if (ip < 0 || ip >= particle->nlocal) {
    if (comm->me == 0) error->message(FLERR, "FixEM: apply_lorentz_force invalid particle index");
    return;
  }
  if (m <= 0.0) {
    if (comm->me == 0) error->message(FLERR, "FixEM: apply_lorentz_force zero or negative mass");
    return;
  }
  if (std::fabs(q) < 1e-12) return;

  Particle::OnePart &p = particle->particles[ip];
  double dt = update->dt;

  if (comm->me == 0 && ip == 0) {
    char msg[512];
    snprintf(msg, 512,
             "FixEM: apply_lorentz_force for particle %d:\n"
             "  pos=(%g,%g,%g)\n"
             "  vel_before=(%g,%g,%g)\n"
             "  E=(%g,%g,%g)\n"
             "  B=(%g,%g,%g)\n"
             "  q=%g, m=%g, dt=%g",
             ip, p.x[0], p.x[1], p.x[2], p.v[0], p.v[1], p.v[2],
             evec[0], evec[1], evec[2], bvec[0], bvec[1], bvec[2], q, m, dt);
    error->message(FLERR, msg);
  }

  double vxb[3] = {
    p.v[1] * bvec[2] - p.v[2] * bvec[1],
    p.v[2] * bvec[0] - p.v[0] * bvec[2],
    p.v[0] * bvec[1] - p.v[1] * bvec[0]
  };

  if (comm->me == 0 && ip == 0) {
    char msg[256];
    snprintf(msg, 256, "FixEM: v x B = (%g, %g, %g)", vxb[0], vxb[1], vxb[2]);
    error->message(FLERR, msg);
  }

  double acc[3];
  for (int d = 0; d < 3; d++) {
    acc[d] = q * (evec[d] + vxb[d]) / m;
  }

  if (comm->me == 0 && ip == 0) {
    char msg[256];
    snprintf(msg, 256, "FixEM: acceleration = (%g, %g, %g)", acc[0], acc[1], acc[2]);
    error->message(FLERR, msg);
  }

  for (int d = 0; d < 3; d++) {
    p.v[d] += acc[d] * dt;
  }

  if (comm->me == 0 && ip == 0) {
    char msg[256];
    snprintf(msg, 256, "FixEM: velocity after update = (%g, %g, %g)",
             p.v[0], p.v[1], p.v[2]);
    error->message(FLERR, msg);
  }
}

void FixEM::set_electric_field(double ex, double ey, double ez)
{ 
  efield[0] = ex; 
  efield[1] = ey; 
  efield[2] = ez; 
  if (field_type < UNIFORM) field_type = UNIFORM; 
  if (comm->me == 0) {
    char msg[128];
    snprintf(msg, 128, "FixEM: set_electric_field to (%g, %g, %g)",
             ex, ey, ez);
    error->message(FLERR, msg);
  }
}

void FixEM::set_magnetic_field(double bx, double by, double bz)
{ 
  bfield[0] = bx; 
  bfield[1] = by; 
  bfield[2] = bz; 
  if (field_type < UNIFORM) field_type = UNIFORM; 
  if (comm->me == 0) {
    char msg[128];
    snprintf(msg, 128, "FixEM: set_magnetic_field to (%g, %g, %g)",
             bx, by, bz);
    error->message(FLERR, msg);
  }
}

void FixEM::set_electric_field_function(FieldFunction func)
{ 
  efunc = func; 
  if (field_type < CUSTOM) field_type = CUSTOM; 
  if (comm->me == 0) error->message(FLERR, "FixEM: set_electric_field_function called");
}

void FixEM::set_magnetic_field_function(FieldFunction func)
{ 
  bfunc = func; 
  if (field_type < CUSTOM) field_type = CUSTOM; 
  if (comm->me == 0) error->message(FLERR, "FixEM: set_magnetic_field_function called");
}

void FixEM::update_custom(int ip, double x, double y, double z, double *fvec)
{
  if (field_type != GRID) return;

  if (ip < 0 || ip >= particle->nlocal) {
    if (comm->me == 0) error->message(FLERR, "FixEM: update_custom invalid particle index");
    return;
  }

  int ic = particle->particles[ip].icell;
  if (ic >= 0 && ic < nglocal) {
    egrid[ic][0] = fvec[0]; 
    egrid[ic][1] = fvec[1]; 
    egrid[ic][2] = fvec[2];
    bgrid[ic][0] = fvec[3]; 
    bgrid[ic][1] = fvec[4]; 
    bgrid[ic][2] = fvec[5];
    if (comm->me == 0) {
      char msg[256];
      snprintf(msg, 256, "FixEM: update_custom updated grid cell %d for particle %d", ic, ip);
      error->message(FLERR, msg);
    }
  } else {
    if (comm->me == 0) {
      char msg[256];
      snprintf(msg, 256, "FixEM: update_custom invalid cell %d for particle %d", ic, ip);
      error->message(FLERR, msg);
    }
  }
}

void FixEM::init_physics()
{ 
  e_charge = 1.602176634e-19;  // Elementary charge in Coulombs
  k_boltz = 1.380649e-23;     // Boltzmann constant in J/K
  mvv2e = 1.0;                // Mass-velocity conversion factor
  if (comm->me == 0) error->message(FLERR, "FixEM: init_physics called");
}

void FixEM::check_species()
{
  species_with_charge = false;
  for (int i = 0; i < particle->nspecies; i++) {
    if (std::fabs(particle->species[i].charge) > 1e-12) { 
      species_with_charge = true; 
      if (particle->species[i].mass <= 0.0) {
        char msg[128];
        snprintf(msg, 128, "FixEM: specie %d has non-positive mass: %g", i, particle->species[i].mass);
        error->all(FLERR, msg);
      }
    }
  }
  if (!species_with_charge && comm->me == 0)
    error->warning(FLERR, "FixEM: no charged species");
  else if (comm->me == 0) {
    char msg[128];
    snprintf(msg, 128, "FixEM: charged species present: %s", species_with_charge ? "yes" : "no");
    error->message(FLERR, msg);
  }
}

void FixEM::allocate()
{
  nglocal = grid->nlocal;
  memory->create(egrid, nglocal, 3, "fix/em:egrid");
  memory->create(bgrid, nglocal, 3, "fix/em:bgrid");

  for (int i = 0; i < nglocal; i++) {
    for (int d = 0; d < 3; d++) {
      egrid[i][d] = 0.0;
      bgrid[i][d] = 0.0;
    }
  }
  if (comm->me == 0) {
    char msg[128];
    snprintf(msg, 128, "FixEM: allocated EM grid with %d cells", nglocal);
    error->message(FLERR, msg);
  }
}

void FixEM::compute_em_grid()
{
  for (int i = 0; i < nglocal; i++) {
    for (int d = 0; d < 3; d++) { 
      egrid[i][d] = efield[d]; 
      bgrid[i][d] = bfield[d]; 
    }
  }
  if (comm->me == 0) error->message(FLERR, "FixEM: compute_em_grid completed");
}

void FixEM::read_efield(const char *filename)
{ 
  if (comm->me == 0) 
    error->warning(FLERR, "FixEM: read_efield not implemented"); 
}

void FixEM::read_bfield(const char *filename)
{ 
  if (comm->me == 0) 
    error->warning(FLERR, "FixEM: read_bfield not implemented"); 
}
