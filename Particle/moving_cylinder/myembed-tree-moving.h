#ifndef BASILISK_HEADER_21
#define BASILISK_HEADER_21
#line 1 "/home/sagar/basilisk/basilisk/src/../myembed-tree-moving.h"
/**
# Embedded boundaries on adaptive trees

This file defines the restriction/prolongation functions which are
necessary to implement [embedded boundaries](embed.h) on adaptive
meshes.

## Volume fraction field *cs*

For the embedded fraction field *cs*, the function below is modelled
closely on the volume fraction refinement function
[fraction_refine()](fractions.h#fraction_refine). */

static void embed_fraction_refine (Point point, scalar cs)
{
  double cc = cs[];

  /**
  If the cell is empty or full, simple injection from the coarse cell
  value is used. */
  
  if (cc <= 0. || cc >= 1.) {
    foreach_child()
      cs[] = cc;
  }
  else {

    /**
    If the cell contains the embedded boundary, we reconstruct the
    boundary using VOF linear reconstruction and a normal estimated
    from the surface fractions. */

    coord n = facet_normal (point, cs, fs);
    double alpha = plane_alpha (cc, n);
      
    foreach_child() {
      static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
      coord nc;
      foreach_dimension()
	nc.x = child.x*n.x;
      cs[] = rectangle_fraction (nc, alpha, a, b);
    }
  }
}

/**
## Surface fractions field *fs*

The embedded surface fractions *fs* are reconstructed using this
function. */

foreach_dimension()
static void embed_face_fraction_refine_x (Point point, scalar s)
{
  vector fs = s.v;

  /**
  If the cell is empty or full, simple injection from the coarse cell
  value is used. */
  
  if (cs[] <= 0. || cs[] >= 1.) {

    /**
    We need to make sure that the fine cells face fractions match
    those of their neighbours. */

    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	fine(fs.x,1,j,k) = cs[];
    for (int i = 0; i <= 1; i++)
      if (!is_refined(neighbor(2*i-1)) && neighbor(2*i-1).neighbors &&
	  (is_local(cell) || is_local(neighbor(2*i-1))))
	for (int j = 0; j <= 1; j++)
	  for (int k = 0; k <= 1; k++)
	    fine(fs.x,2*i,j,k) = fs.x[i];
  }
  else {

    /**
    If the cell contains the embedded boundary, we reconstruct the
    boundary using VOF linear reconstruction and a normal estimated
    from the surface fractions. */

    coord n = facet_normal (point, cs, fs);
    double alpha = plane_alpha (cs[], n);
      
    /**
    We need to reconstruct the face fractions *fs* for the fine cells.
    
    For the fine face fractions contained within the coarse cell,
    we compute the intersections directly using the VOF
    reconstruction. */

#if dimension == 2

    /**
    In 2D, we obtain the face fractions by taking into
    account the orientation of the normal. */

    if (2.*fabs(alpha) < fabs(n.y)) {
      double yc = alpha/n.y;
      int i = yc > 0.;
      fine(fs.x,1,1 - i) = n.y < 0. ? 1. - i : i;
      fine(fs.x,1,i) = n.y < 0. ? i - 2.*yc : 1. - i + 2.*yc;
    }
    else
      fine(fs.x,1,0) = fine(fs.x,1,1) = alpha > 0.;

#else // dimension == 3

    /**
    in 3D, we use the 2D projection of the reconstruction. */

    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	// fixme: change to avoid accesing *cs* in non-initialized cells
	// when running on multiple procs (-catch)
	if (!fine(cs,0,j,k) || !fine(cs,1,j,k))
	  fine(fs.x,1,j,k) = 0.;
	else {
	  static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
	  coord nc;
	  nc.x = 0., nc.y = (2.*j - 1.)*n.y, nc.z = (2.*k - 1.)*n.z;
	  fine(fs.x,1,j,k) = rectangle_fraction (nc, alpha, a, b);
	}

#endif // dimension == 3
    
    /**
    For the fine face fractions coincident with the faces of the
    coarse cell, we obtain the intersection position from the
    coarse cell face fraction. */

    for (int i = 0; i <= 1; i++)
      if (neighbor(2*i-1).neighbors &&
	  (is_local(cell) || is_local(neighbor(2*i-1)))) {
	if (!is_refined(neighbor(2*i-1))) {
	  if (fs.x[i] <= 0. || fs.x[i] >= 1.)
	    for (int j = 0; j <= 1; j++)
	      for (int k = 0; k <= 1; k++)
		fine(fs.x,2*i,j,k) = fs.x[i];
	  else {
#if dimension == 2
	  
	    /**
	    In 2D the orientation is obtained by looking at the values
	    of face fractions in the transverse direction. */
	  
	    double a = fs.y[0,1] <= 0. || fs.y[2*i-1,1] <= 0. ||
	      fs.y[] >= 1. || fs.y[2*i-1] >= 1.;
	    if ((2.*a - 1)*(fs.x[i] - 0.5) > 0.) {
	      fine(fs.x,2*i,0) = a;
	      fine(fs.x,2*i,1) = 2.*fs.x[i] - a;
	    }
	    else {
	      fine(fs.x,2*i,0) = 2.*fs.x[i] + a - 1.;
	      fine(fs.x,2*i,1) = 1. - a;
	    }

#else  // dimension == 3

	    /**
	    In 3D we reconstruct the face fraction from the projection
	    of the cell interface reconstruction, as above. */
	  
	    for (int j = 0; j <= 1; j++)
	      for (int k = 0; k <= 1; k++) {
		static const coord a = {0.,0.,0.}, b = {.5,.5,.5};
		coord nc;
		nc.x = 0., nc.y = (2.*j - 1.)*n.y, nc.z = (2.*k - 1.)*n.z;
		fine(fs.x,2*i,j,k) =
		  rectangle_fraction (nc, alpha - n.x*(2.*i - 1.)/2., a, b);
	      }

#endif // dimension == 3
	  }
	}

	/**
	The face fractions of empty children cells must be zero. */
	
	for (int j = 0; j <= 1; j++)
	#if dimension > 2
	  for (int k = 0; k <= 1; k++)
	#endif
	    if (fine(fs.x,2*i,j,k) && !fine(cs,i,j,k))
	      fine(fs.x,2*i,j,k) = 0.;
      }
  }
}

/**
## Restriction of cell-centered fields

We now define the restriction function for cell-centered fields. The
goal is to define second-order operators which do not use any values
from cells entirely contained within the embedded boundary (for which
*cs = 0*) as well as, if possible, values from emerged cells (*csm1[]
= 0*) that have not yet been updated. Without loss of generality, we
consider a quadtree and proceed as follows:

* if four children cells are available, we use the simple diagonal
averaging:
$$
s\left[\,\right] =
\frac{1}{2}\left(\frac{\mathrm{fine}\left(s\left[\,\right]\right) +
\mathrm{fine}\left(s\left[1,1\right]\right)}{2} +
\frac{\mathrm{fine}\left(s\left[1\right]\right) +
\mathrm{fine}\left(s\left[0,1\right]\right)}{2}\right);
$$

* if three children cells are available, we also use diagonal
averaging, but only along the available diagonal:
$$
s\left[\,\right] = \left\{
\begin{aligned}
&
\mathrm{if}\: \mathrm{fine}\left(c\left[\,\right]\right) > 0 \:
\mathrm{and} \: \mathrm{fine}\left(c \left[1,1\right]\right) > 0 :\\
& 
\qquad \frac{\mathrm{fine}\left(s\left[\,\right]\right) +
\mathrm{fine}\left(s\left[1,1\right]\right)}{2} \\    
& 
\mathrm{if}\: \mathrm{fine}\left(c\left[1\right]\right) > 0 \:
\mathrm{and} \: \mathrm{fine}\left(c \left[0,1\right]\right) > 0 :\\
&
\qquad \frac{\mathrm{fine}\left(s\left[1\right]\right) +
\mathrm{fine}\left(s\left[0,1\right]\right)}{2}.
\end{aligned}
\right.
$$

When restricting it is unfortunately not always possible to obtain a
second-order interpolation. If less than three children cells are
available in 2D, we face such a pathological situation. In this case,
some external information (i.e. a boundary gradient condition) is
required to be able to maintain second-order accuracy. This
information can be passed by defining the *embed_gradient()* function
of the field being restricted. In this case, we first compute the
average value of $s$, denoted $\bar{s}$, over all available children
cells as well as the barycenter $\mathbf{b}$ of the available children
cells. We then use, if possible, the user-provided value of the
gradient of $s$ defined at the center of the cell, denoted
$\mathbf{\nabla} s\left[\,\right]$, and a Taylor expansion to improve
our approximation of $s$:
$$
s[] = \bar{s} - \mathbf{b}\cdot\mathbf{\nabla} s\left[\,\right].
$$
If an a priori value of $\mathbf{\nabla} s\left[\,\right]$ is
not available, the restriction function in this pathological situation
is first-order only.

![Graphical representation of the restriction of a cell-centered
scalar on a 2D grid.](basilisk-embed-restriction.jpeg){width=100%} */

attribute {
  void (* embed_gradient) (Point, scalar, coord *);
}

static inline void restriction_embed_linear (Point point, scalar s)
{  
  // 0 children
  if (!cs[]) {
    s[] = 0.;
    return;
  }

  /**
  We first try to interpolate "diagonally". If enough child cells are
  defined (i.e. have non-zero embedded fractions and are not emerged
  cells that have not yet been updated), we return the corresponding
  value. */

  double val = 0., nv = 0.;
  for (int i = 0; i <= 1; i++)
#if dimension > 2
    for (int j = 0; j <= 1; j++)
#endif
      if (fine(cs,0,i,j) && fine(cs,1,!i,!j) &&
	  (emerged || (fine(csm1,0,i,j) && fine(csm1,1,!i,!j))))
	val += (fine(s,0,i,j) + fine(s,1,!i,!j))/2., nv++;
  if (nv > 0.) {
    s[] = val/nv;
    return;
  }

  /**
  Otherwise, we use the average of the child cells which are
  defined. In the case of a fixed embedded boundary, there is at least
  one. However, for moving boundaries, if all child cells are emerged
  cells, we may use values which are not properly defined. */
  
  coord p = {0.,0.,0.};
  foreach_child()
    if (cs[] && (emerged || csm1[]))
      p.x += x, p.y += y, p.z += z, val += s[], nv++;
  if (nv > 0.)
    s[] = val/nv;
  else { //We resort to using emerged cell values
    foreach_child()
      if (cs[])
	p.x += x, p.y += y, p.z += z, val += s[], nv++;
    assert (nv > 0.);
    s[] = val/nv;
  }

  /**
  If the gradient is defined and if the variable is not using
  homogeneous boundary conditions (i.e. the variable is not a
  "correction" used by the Poisson solver), we improve the
  interpolation using this information. */

  bool homogeneous = true;
  for (int b = 0; b < nboundary; b++)
    if (s.boundary[b] != s.boundary_homogeneous[b])
      homogeneous = false;

  if (s.embed_gradient && !homogeneous) {
    coord o = {x,y,z}, g;
    s.embed_gradient (point, s, &g);
    foreach_dimension()
      s[] += (o.x - p.x/nv)*g.x;
  }
}

/**
## Refinement/prolongation of cell-centered fields

For refinement, we use either bilinear interpolation, if the required
four coarse cell values are defined or trilinear interpolation if only
three coarse cell values are defined. If less that three coarse cell
values are defined ("pathological cases" below), we try to estimate
gradients in each direction and add the corresponding
correction. Without loss of generality, we consider a quadtree and
proceed as follows:

* if all four coarse cells are available, we use the bilinear
interpolation:
$$
s\left[\,\right] = \frac{9 \, \mathrm{coarse}
\left(s\left[\,\right]\right) +
3\left(\mathrm{coarse}\left(s\left[\mathrm{child}.x\right]\right) +
\mathrm{coarse}\left(s\left[0,\mathrm{child}.y\right]\right)\right) +
\mathrm{coarse}\left(s\left[\mathrm{child}.x,\mathrm{child}.y\right]\right)}{16},
$$

* if three coarse cells are available, excluding the diagonal coarse
cell $\left[\mathrm{child}.x,\mathrm{child}.y\right]$, we use the
following triangular interpolation:
$$
s\left[\,\right] =  \frac{2 \, \mathrm{coarse}
\left(s\left[\,\right]\right) + \mathrm{coarse}\left(s\left[\mathrm{child}.x\right]\right) +
\mathrm{coarse}\left(s\left[0,\mathrm{child}.y\right]\right)}{4}.
$$

* if only three coarse cells are available, including the diagonal
coarse cell $\left[\mathrm{child}.x,\mathrm{child}.y\right]$, we use
the following diagonal interpolation:
$$
s\left[\,\right] = \frac{3 \, \mathrm{coarse}
\left(s\left[\,\right]\right) +
\mathrm{coarse}\left(s\left[\mathrm{child}.x,\mathrm{child}.y\right]\right)}{4}.
$$

* if only one or two coarse cells are available, we face a
pathological situation. In this case, we compute the value of $s$
using a Taylor expansion, where, if possible, one simple face gradient
$\bar{\nabla}_{d}^{f} s$ at level $l$ is computed per dimension,
preferably along the faces of the parent cell that follow the natural
orientation of the halo cell with respect to its parent cell. For
example:
$$
s\left[\,\right] = \mathrm{coarse} \left(s\left[\,\right]\right) +
\frac{\Delta}{4}\mathrm{coarse} \left(\bar{\nabla}_{x}^{f} s
\left[1\right]\right) + \frac{\Delta}{4} \mathrm{coarse}
\left(\bar{\nabla}_{y}^{f} s \left[\,\right]\right).
$$

![Graphical representation of the prolongation of a cell-centered
scalar on a 2D grid.](basilisk-embed-prolongation.jpeg){width=100%} */

static inline void refine_embed_linear (Point point, scalar s)
{  
  foreach_child() {
    if (!cs[])
      s[] = 0.;
    else {
      assert (coarse(cs));
      int i = (child.x + 1)/2, j = (child.y + 1)/2;
#if dimension == 2
      if (coarse(fs.x,i) > 0.25 && coarse(fs.y,0,j) > 0.25 &&
	  (coarse(cs) == 1. || coarse(cs,child.x) == 1. ||
	   coarse(cs,0,child.y) == 1. || coarse(cs,child.x,child.y) == 1.) &&
	  (emerged || (coarse(csm1) && coarse(csm1,child.x) &&
		       coarse(csm1,0,child.y) && coarse(csm1,child.x,child.y)))) {
	assert (coarse(cs,child.x) && coarse(cs,0,child.y));
	if (coarse(fs.x,i,child.y) && coarse(fs.y,child.x,j)) {
	  // bilinear interpolation
	  assert (coarse(cs,child.x,child.y));
	  s[] = (9.*coarse(s) + 
		 3.*(coarse(s,child.x) + coarse(s,0,child.y)) + 
		 coarse(s,child.x,child.y))/16.;
	}
	else
	  // triangular interpolation	  
	  s[] = (2.*coarse(s) + coarse(s,child.x) + coarse(s,0,child.y))/4.;
      }
      else if (coarse(cs,child.x,child.y) &&
	       ((coarse(fs.x,i) && coarse(fs.y,child.x,j)) ||
		(coarse(fs.y,0,j) && coarse(fs.x,i,child.y))) &&
	       (emerged || (coarse(csm1) && coarse(csm1,child.x,child.y)))) {
	// diagonal interpolation
	s[] = (3.*coarse(s) + coarse(s,child.x,child.y))/4.;
      }
#else // dimension == 3
      int k = (child.z + 1)/2;
      if (coarse(fs.x,i) > 0.25 && coarse(fs.y,0,j) > 0.25 &&
	  coarse(fs.z,0,0,k) > 0.25 &&
	  (coarse(cs) == 1. || coarse(cs,child.x) == 1. ||
	   coarse(cs,0,child.y) == 1. || coarse(cs,child.x,child.y) == 1. ||
	   coarse(cs,0,0,child.z) == 1. || coarse(cs,child.x,0,child.z) == 1. ||
	   coarse(cs,0,child.y,child.z) == 1. ||
	   coarse(cs,child.x,child.y,child.z) == 1.) &&
	  (emerged || (coarse(csm1) && coarse(csm1,child.x) &&
		       coarse(csm1,0,child.y) && coarse(csm1,child.x,child.y) &&
		       coarse(csm1,0,0,child.z) && coarse(csm1,child.x,0,child.z) &&
		       coarse(csm1,0,child.y,child.z) && coarse(csm1,child.x,child.y,child.z)))) {
	assert (coarse(cs,child.x) && coarse(cs,0,child.y) &&
		coarse(cs,0,0,child.z));
	if (coarse(fs.x,i,child.y) && coarse(fs.y,child.x,j) &&
	    coarse(fs.x,i,0,child.z) && coarse(fs.y,0,j,child.z) &&
	    coarse(fs.x,i,child.y,child.z) && coarse(fs.y,child.x,j,child.z) &&
	    //
	    coarse(fs.z,child.x,child.y,k) &&
	    coarse(fs.z,child.x,0,k) && coarse(fs.z,0,child.y,k)) {
	  assert (coarse(cs,child.x,child.y) && coarse(cs,child.x,0,child.z) &&
		  coarse(cs,0,child.y,child.z) &&
		  coarse(cs,child.x,child.y,child.z));
	  // bilinear interpolation
	  s[] = (27.*coarse(s) + 
		 9.*(coarse(s,child.x) + coarse(s,0,child.y) +
		     coarse(s,0,0,child.z)) + 
		 3.*(coarse(s,child.x,child.y) + coarse(s,child.x,0,child.z) +
		     coarse(s,0,child.y,child.z)) + 
		 coarse(s,child.x,child.y,child.z))/64.;
	}
	else
	  // tetrahedral interpolation
	  s[] = (coarse(s) + coarse(s,child.x) + coarse(s,0,child.y) +
	  	 coarse(s,0,0,child.z))/4.;
      }
      else if (coarse(cs,child.x,child.y,child.z) &&
	       ((coarse(fs.z,child.x,child.y,k) &&
		 ((coarse(fs.x,i) && coarse(fs.y,child.x,j)) ||
		  (coarse(fs.y,0,j) && coarse(fs.x,i,child.y))))
		||
		(coarse(fs.z,0,0,k) &&
		 ((coarse(fs.x,i,0,child.z) && coarse(fs.y,child.x,j,child.z)) ||
		  (coarse(fs.y,0,j,child.z) && coarse(fs.x,i,child.y,child.z))))
		||
		(coarse(fs.z,child.x,0,k) &&
		 coarse(fs.x,i) && coarse(fs.y,child.x,j,child.z))
		||
		(coarse(fs.z,0,child.y,k) &&
		 coarse(fs.y,0,j) && coarse(fs.x,i,child.y,child.z))
		) &&
	       (emerged || (coarse(csm1) && coarse(csm1,child.x,child.y,child.z))))
	// diagonal interpolation
	s[] = (3.*coarse(s) + coarse(s,child.x,child.y,child.z))/4.;
#endif // dimension == 3
      else {
	// Pathological cases, use 1D gradients.
	s[] = coarse(s);
	foreach_dimension() {
	  if (coarse(fs.x,(child.x + 1)/2) && coarse(cs,child.x) &&
	      (emerged || (coarse(csm1) && coarse(csm1,child.x))))
	    s[] += (coarse(s,child.x) - coarse(s))/4.;
	  else if (coarse(fs.x,(- child.x + 1)/2) && coarse(cs,- child.x) &&
		   (emerged || (coarse(csm1) && coarse(csm1,- child.x))))
	    s[] -= (coarse(s,- child.x) - coarse(s))/4.;
	}
      }
    }
  }
}

/**
## Refinement/prolongation of face-centered velocity

This function is modelled on
[*refine_face_x()*](/src/grid/tree-common.h#refine_face_x) and is
typically used to refine the values of the face-centered velocity
field *uf*. It uses linear interpolation, taking into account the
weighting by the embedded fractions *fs*. */

foreach_dimension()
void refine_embed_face_x (Point point, scalar s)
{
  vector v = s.v;
  for (int i = 0; i <= 1; i++)
    if (neighbor(2*i - 1).neighbors &&
	(is_local(cell) || is_local(neighbor(2*i - 1)))) {
      double g1 = fs.x[i] >= 1. && fs.x[i,+1] && fs.x[i,-1] ?
	(v.x[i,+1]/fs.x[i,+1] - v.x[i,-1]/fs.x[i,-1])/8. : 0.;
      double g2 = fs.x[i] >= 1. && fs.x[i,0,+1] && fs.x[i,0,-1] ?
	(v.x[i,0,+1]/fs.x[i,0,+1] - v.x[i,0,-1]/fs.x[i,0,-1])/8. : 0.;
      for (int j = 0; j <= 1; j++)
	for (int k = 0; k <= 1; k++)
	  fine(v.x,2*i,j,k) = fs.x[i] ?
	    fine(fs.x,2*i,j,k)*(v.x[i]/fs.x[i] +
				(2*j - 1)*g1 + (2*k - 1)*g2) : 0.;
    }
  if (is_local(cell)) {
    double g1 = (fs.x[0,+1] + fs.x[1,+1]) && (fs.x[0,-1] + fs.x[1,-1]) ?
      ((v.x[0,+1] + v.x[1,+1])/(fs.x[0,+1] + fs.x[1,+1]) -
       (v.x[0,-1] + v.x[1,-1])/(fs.x[0,-1] + fs.x[1,-1]))/8. : 0.;
    double g2 = (fs.x[1,0,+1] + fs.x[0,0,+1]) && (fs.x[1,0,-1] + fs.x[0,0,-1]) ?
      ((v.x[0,0,+1] + v.x[1,0,+1])/(fs.x[1,0,+1] + fs.x[0,0,+1]) -
       (v.x[0,0,-1] + v.x[1,0,-1])/(fs.x[1,0,-1] + fs.x[0,0,-1]))/8. : 0.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
	fine(v.x,1,j,k) = fs.x[] + fs.x[1] ?
	  fine(fs.x,1,j,k)*((v.x[] + v.x[1])/(fs.x[] + fs.x[1]) +
			    (2*j - 1)*g1 + (2*k - 1)*g2) : 0.;
  }
}

#endif
