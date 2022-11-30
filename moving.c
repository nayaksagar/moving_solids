#include "navier-stokes/centered.h"
#include "fractions.h"
//#include "embed.h"
// face vector muc[];

int main(){
  L0 = 8.;
  origin (-0.5, -L0/2.);
  N = 512;
  const face vector muc[] = {0.00078125,0.00078125};
  mu = muc;
  run();
}

// event properties(i++){
//   foreach_face()
//     muv.x[]=fm.x[]*0.125/Re;
// }

event init (t = 0) {
  mask(y > 0.5 ? top: y < -0.5 ? bottom : none);
}

scalar cylinder[];

event moving_cylinder (i++) {
  coord vc = {1.*t,0.}; // the velocity of the cylinder
  fraction (cylinder, - (sq(x - vc.x*t) + sq(y - vc.y*t) - sq(0.0625)));


  foreach()
    foreach_dimension()
      u.x[] = cylinder[]*vc.x + (1. - cylinder[])*u.x[];
  boundary ((scalar *){u});
}

event logfile (i++)
  fprintf (stderr, "%d %g %d %d %d %d\n", i, t,
	   mgp.i, mgp.nrelax, mgu.i, mgu.nrelax);

event images (t += 0.1; t <= 10.0) {
  static FILE * fp = popen ("ppm2gif > vort.gif", "w");
  scalar omega[];
  vorticity (u, omega);

  scalar m[];
  foreach()
    m[] = 0.5 - cylinder[];
  boundary ({m});

  output_ppm (omega, fp, box = {{-0.5,-0.5},{7.5,0.5}}, mask = m,
	      min=-10, max=10, linear=true);
}

event adapt (i++) {
  adapt_wavelet ((scalar *){u}, (double[]){3e-2,3e-2}, 9, 4);
}// event logfile(i+=10){
//   coord Fp, Fmu;
//   embed_force (p, u, mu, &Fp, &Fmu);
//   fprintf (stderr, "%d %g %g %.6f %.6f %.6f %.6f\n",
// 	   i, t, dt, Fp.x, Fp.y, Fmu.x, Fmu.y);
// }

// event logfile(i+=10){
//   coord Fp, Fmu;
//   embed_force (p, u, mu, &Fp, &Fmu);
//   fprintf (stderr, "%d %g %g %.6f %.6f %.6f %.6f\n",
// 	   i, t, dt, Fp.x, Fp.y, Fmu.x, Fmu.y);
// }
