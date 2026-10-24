"""Exercise actual cooling dust functions with controlled unit/hydro adapters."""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

root = Path(__file__).resolve().parents[1]
source = (root / "src/cooling/KIARA/cooling.c").read_text()


def extract(name):
    prefix = {"cooling_sputter_dust": "__attribute__((always_inline)) INLINE void ",
              "cooling_init_chemistry": "void "}.get(name, "static void ")
    signature = source.index(prefix + name + "(")
    brace = source.index("{", signature)
    depth = 1
    end = brace + 1
    while depth:
        depth += (source[end] == "{") - (source[end] == "}")
        end += 1
    return source[signature:end]


prefix = r'''
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "chemistry/KIARA/dust_sources.h"
#include "chemistry/KIARA/dust_evolution.h"
#define INLINE static inline
#define error(...) do { fprintf(stderr, __VA_ARGS__); abort(); } while (0)
#define warning(...) ((void)0)
enum {chemistry_element_H, chemistry_element_He, chemistry_element_C,
      chemistry_element_N, chemistry_element_O, chemistry_element_Ne,
      chemistry_element_Mg, chemistry_element_Si, chemistry_element_S,
      chemistry_element_Ca, chemistry_element_Fe, chemistry_element_count};
enum {UNIT_CONV_TIME, UNIT_CONV_DENSITY, UNIT_CONV_TEMPERATURE, UNIT_CONV_MASS};
struct unit_system {double factor[4];};
struct cosmology {double a;};
struct cooling_function_data {
  int use_grackle_dust_evol;
  double dust_grain_sizes[KIARA_DUST_N_BINS], dust_initial_distribution[KIARA_DUST_N_BINS],
         dust_growth_tauref, dust_growth_densref, dust_growth_Tref,
         dust_sne_shockspeed, dust_sne_coeff, dust_destruction_eff;
  double self_enrichment_metallicity;
  struct {double SolarAbundances[10];} chemistry;
};
struct cooling_part_data {
  float dust_mass, subgrid_temp, dust_mass_fraction[11],
        dust_size_distribution[KIARA_DUST_N_MATERIALS][KIARA_DUST_N_BINS];
};
struct part {
  long long id;
  float mass, rho, u;
  struct cooling_part_data cooling_data;
  struct {float metal_mass_fraction[11], metal_mass_fraction_total;} chemistry_data;
  struct {float SNe_ThisTimeStep;} feedback_data;
};
struct xpart {struct {double e_frac;} cooling_data;};
static double units_cgs_conversion_factor(const struct unit_system *u, int key) {return u->factor[key];}
static double hydro_get_mass(const struct part *p) {return p->mass;}
static double hydro_get_physical_density(const struct part *p, const struct cosmology *c) {return p->rho / pow(c->a, 3.);}
static double hydro_get_physical_internal_energy(const struct part *p, const struct xpart *x, const struct cosmology *c) {(void)x; return p->u / (c->a*c->a);}
static double cooling_convert_u_to_temp(double u, double e, const struct cooling_function_data *c, const struct part *p, const struct xpart *x) {(void)e;(void)c;(void)p;(void)x;return u;}
static double chemistry_get_total_metal_mass_fraction_for_cooling(const struct part *p) {return p->chemistry_data.metal_mass_fraction_total;}
'''
body = r'''
static struct part particle(const struct cooling_function_data *cooling,
    const struct unit_system *us, double a, double temperature) {
  struct part p={0}; p.mass=1.f;
  p.rho=1.e-27*pow(a,3.)/us->factor[UNIT_CONV_DENSITY];
  p.u=temperature*a*a/us->factor[UNIT_CONV_TEMPERATURE];
  p.cooling_data.dust_mass=.01f;
  p.cooling_data.dust_mass_fraction[0]=1.f;
  p.cooling_data.dust_mass_fraction[2]=.3f;
  p.cooling_data.dust_mass_fraction[4]=.7f;
  p.cooling_data.dust_size_distribution[0][0]=.3f;
  p.cooling_data.dust_size_distribution[1][KIARA_DUST_N_BINS-1]=.7f;
  p.chemistry_data.metal_mass_fraction_total=.01f;
  p.chemistry_data.metal_mass_fraction[0]=.74f;
  p.chemistry_data.metal_mass_fraction[1]=.24f;
  p.chemistry_data.metal_mass_fraction[2]=.005f;
  p.chemistry_data.metal_mass_fraction[4]=.005f;
  (void)cooling;
  return p;
}
static double inventory(const struct part *p, int elem) {
  return p->mass*p->chemistry_data.metal_mass_fraction[elem] +
         p->cooling_data.dust_mass*p->cooling_data.dust_mass_fraction[elem];
}
static void check_state(const struct part *p) {
  double bins=0., composition=0.;
  for (int mat=0;mat<2;++mat)
    for(int i=0;i<KIARA_DUST_N_BINS;++i) {
      assert(p->cooling_data.dust_size_distribution[mat][i]>=0.f);
      bins+=p->cooling_data.dust_size_distribution[mat][i];
    }
  for(int e=1;e<11;++e) composition+=p->cooling_data.dust_mass_fraction[e];
  assert(fabs(bins-(p->cooling_data.dust_mass>0.f))<1.e-6);
  assert(fabs(composition-(p->cooling_data.dust_mass>0.f))<1.e-6);
  assert(fabs(p->chemistry_data.metal_mass_fraction[0]+
      p->chemistry_data.metal_mass_fraction[1]+
      p->chemistry_data.metal_mass_fraction_total+
      p->cooling_data.dust_mass/p->mass-1.)<2.e-7);
  double carbon=0.;
  for(int i=0;i<KIARA_DUST_N_BINS;++i)
    carbon+=p->cooling_data.dust_size_distribution[0][i];
  assert(fabs(carbon-p->cooling_data.dust_mass_fraction[2])<1.e-6);
}
int main(void) {
  struct cooling_function_data cooling={.use_grackle_dust_evol=1,
      .dust_growth_tauref=1.,
      .dust_growth_densref=1.e-27,.dust_growth_Tref=20.,
      .dust_sne_shockspeed=100.,.dust_sne_coeff=1.,.dust_destruction_eff=.3};
  assert(kiara_dust_grid_logarithmic(1.e-4,1.,cooling.dust_grain_sizes));
  assert(kiara_dust_source_lognormal(cooling.dust_grain_sizes,.1,.6,0.,cooling.dust_initial_distribution));
  struct xpart xp={0};
  for(int epoch=0;epoch<3;++epoch) for(int unit=0;unit<2;++unit) {
    struct cosmology cosmo={.a=(double[]){1.,.5,.1}[epoch]};
    struct unit_system us={{unit ? 3.15569251e16 : 1.,1.e-30,unit ? 1000. : 1.,1.989e33}};
    struct part p=particle(&cooling,&us,cosmo.a,1.e7);
    const double old_c=inventory(&p,2),old_o=inventory(&p,4);
    const double delta=.5*cooling.dust_grain_sizes[0];
    const double tau=1.7e8*3.15569251e7*(pow(2.e6/1.e7,2.5)+1.)/us.factor[0];
    cooling_sputter_dust(&us,&cosmo,&cooling,&p,&xp,delta/.1*tau);
    const double expected=.01*(.3*.125+.7*pow(1.-delta,3.));
    assert(fabs(p.cooling_data.dust_mass/expected-1.)<2.e-6);
    assert(fabs(inventory(&p,2)-old_c)<2.e-9);
    assert(fabs(inventory(&p,4)-old_o)<2.e-9);
    check_state(&p);
    p=particle(&cooling,&us,cosmo.a,1.e7);
    cooling_sputter_dust(&us,&cosmo,&cooling,&p,&xp,2./.1*tau);
    assert(p.cooling_data.dust_mass==0.f);
    assert(fabs(inventory(&p,2)-old_c)<2.e-9);
    assert(fabs(inventory(&p,4)-old_o)<2.e-9);
    check_state(&p);
    p=particle(&cooling,&us,cosmo.a,20.);
    p.cooling_data.dust_size_distribution[1][KIARA_DUST_N_BINS-1]=0.f;
    p.cooling_data.dust_size_distribution[1][0]=.7f;
    const double growth_tau=3.15569251e16/us.factor[0]*(.0134/.01);
    cooling_evolve_grain_size(&us,&cosmo,&cooling,&p,&xp,.03/.1*growth_tau);
    assert(p.cooling_data.dust_mass==.01f);
    double number=0.;
    for(int m=0;m<2;++m) for(int i=0;i<KIARA_DUST_N_BINS;++i)
      number+=p.cooling_data.dust_size_distribution[m][i]/pow(cooling.dust_grain_sizes[i],3.);
    assert(fabs(number*pow(.0301,3.)-1.)<2.e-6);
    check_state(&p);
    p=particle(&cooling,&us,cosmo.a,20.);
    kiara_dust_distribution_zero(p.cooling_data.dust_size_distribution);
    cooling_init_dust_distribution(&cooling,&p);
    for(int i=0;i<KIARA_DUST_N_BINS;++i) {
      assert(fabs(p.cooling_data.dust_size_distribution[0][i]-.3*cooling.dust_initial_distribution[i])<1.e-7);
      assert(fabs(p.cooling_data.dust_size_distribution[1][i]-.7*cooling.dust_initial_distribution[i])<1.e-7);
    }
    check_state(&p);
    p=particle(&cooling,&us,cosmo.a,20.);
    p.chemistry_data.metal_mass_fraction_total=0.f;
    kiara_dust_distribution_zero(p.cooling_data.dust_size_distribution);
    p.cooling_data.dust_size_distribution[0][KIARA_DUST_N_BINS-1]=.3f;
    p.cooling_data.dust_size_distribution[1][KIARA_DUST_N_BINS-1]=.7f;
    p.feedback_data.SNe_ThisTimeStep=(-3.*log(.9))/(.3*6800.);
    cooling_evolve_grain_size(&us,&cosmo,&cooling,&p,&xp,1.);
    number=0.;
    for(int m=0;m<2;++m) for(int i=0;i<KIARA_DUST_N_BINS;++i)
      number+=p.cooling_data.dust_size_distribution[m][i]/pow(cooling.dust_grain_sizes[i],3.);
    assert(fabs(number*pow(.9,3.)-1.)<2.e-6);
    assert(p.cooling_data.dust_mass==.01f);
  }
  cooling.self_enrichment_metallicity=1.e-7;
  cooling.chemistry.SolarAbundances[0]=.24;
  cooling.chemistry.SolarAbundances[1]=.001;
  cooling.chemistry.SolarAbundances[3]=.004;
  for(int agb=0;agb<2;++agb) {
    assert(kiara_dust_source_lognormal(cooling.dust_grain_sizes,.1,
        agb ? .47 : .6,agb ? 4. : 0.,cooling.dust_initial_distribution));
    struct part seed={0}; seed.mass=1.f; seed.cooling_data.subgrid_temp=20.f;
    cooling_init_chemistry(&cooling,&seed);
    assert(fabs(seed.cooling_data.dust_mass/2.e-8-1.)<2.e-7);
    assert(fabs(seed.chemistry_data.metal_mass_fraction_total/8.e-8-1.)<2.e-7);
    for(int i=0;i<KIARA_DUST_N_BINS;++i) {
      assert(fabs(seed.cooling_data.dust_size_distribution[0][i]-.2*cooling.dust_initial_distribution[i])<1.e-7);
      assert(fabs(seed.cooling_data.dust_size_distribution[1][i]-.8*cooling.dust_initial_distribution[i])<1.e-7);
    }
    check_state(&seed);
  }
  printf("PASS actual cooling dust functions: %d bins\n",KIARA_DUST_N_BINS);
  return 0;
}
'''
program = prefix + "\n".join(extract(name) for name in (
    "cooling_init_dust_distribution", "cooling_sputter_dust",
    "cooling_evolve_grain_size", "cooling_init_chemistry")) + body
with tempfile.TemporaryDirectory(prefix="kiara-dust-cooling-") as scratch:
    cfile = Path(scratch) / "cooling_probe.c"
    cfile.write_text(program)
    for count in (2, 3, 6, 12, 41):
        binary = Path(scratch) / f"cooling-{count}"
        subprocess.run(shlex.split(os.environ.get("CC", "cc")) + [
            "-std=c11", "-Wall", "-Wextra", "-Werror", "-fsanitize=address,undefined",
            f"-DKIARA_DUST_N_BINS={count}", "-I", str(root / "src"),
            str(cfile), "-lm", "-o", str(binary)], check=True)
        subprocess.run([str(binary)], check=True)
