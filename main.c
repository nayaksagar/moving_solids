#include "embed.h"
#include "navier-stokes/centered.h"

int maxlevel=9;
double Re=160.;
face vector muv[];

int main(){
  size (8.0);
  origin(-0.5, -L0/2.);
  N=256;
  mu=muv;
  run();
}

event properties(i++){
  foreach_face()
    muv.x[]=fm.x[]*0.125/Re;
}

u.n[left] = dirichlet(0.5);
p[left] = neumann(0.);
pf[left] = neumann(0.);

u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

u.n[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);

event init(t=0){
  solid (cs, fs, intersection (intersection (0.5 - y, 0.5 + y), sq(x) + sq(y) - sq(0.125/2.)));
  foreach()
     u.x[] = cs[] ? 1. : 0.;
}

// event cylinder(i++){
//   coord vc={1.0,0.};
//   solid (cs, fs, intersection (intersection (0.5 - y, 0.5 + y), sq(x-vc.x*t) + sq(y-vc.y*t) - sq(0.125/2.)));
//   foreach()
//     foreach_dimension()
//       u.x[] = cs[]*vc.x + (1. - cs[])*u.x[];
//   boundary ((scalar *){u});
// }

// event logfile (i++)
//   fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

event output(i+=100){
	output_ppm(u.x,box={{-0.5,-0.5},{7.5,0.5}},linear=true);
	// output_ppm(f,file="f.mp4",linear=false);
}

event movies (i += 4; t <= 15.)
  {
    scalar omega[], m[];
    vorticity (u, omega);
    foreach()
      m[] = cs[] - 0.5;
    output_ppm(omega,box={{-0.5,-0.5},{7.5,0.5}},linear=true);
    // output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
  	//       min = -10, max = 10, linear = true, mask = m);
    // output_ppm (f, file = "f.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
  	//       linear = false, min = 0, max = 1, mask = m);
  }

event adapt (i++) {
    adapt_wavelet ({cs,u}, (double[]){1e-2,3e-2,3e-2}, maxlevel, 4);
  }

event logfile(i+=10){
  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);
  fprintf (stderr, "%d %g %g %.6f %.6f %.6f %.6f\n",
	   i, t, dt, Fp.x, Fp.y, Fmu.x, Fmu.y);
}
