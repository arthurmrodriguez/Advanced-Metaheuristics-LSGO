#include "cec2013_funcs.h"

#include "load_data.h"

// Common global vars to every function. Some functions will use
// all the vars, others only a few.
arma::mat xopt;
long double lb, ub;

arma::mat  R100, R50, R25;
arma::uvec p;
arma::vec  s, w;

arma::vec c;
long double m;

// ///////////// //
// Aux functions //
// ///////////// //

void init_simple(mat_t* matfp) {
  xopt=loadMatrix(matfp, (char*) "xopt");
  //xopt=xopt.submat(0,0,9,0);

  arma::mat tmp;
  tmp=loadMatrix(matfp, (char*) "lb");
  lb=tmp(0, 0);
  tmp=loadMatrix(matfp, (char*) "ub");
  ub=tmp(0, 0);
}

void init_extended(mat_t* matfp) {
  init_simple(matfp);

  R100=loadMatrix(matfp, (char*) "R100");
  R50 =loadMatrix(matfp, (char*) "R50" );
  R25 =loadMatrix(matfp, (char*) "R25" );

  p=loadIndices(matfp, (char*) "p");

  arma::mat tmp;
  tmp=loadMatrix(matfp, (char*) "s");
  arma::vec tmp_s(tmp);
  s=tmp_s;
  tmp=loadMatrix(matfp, (char*) "w");
  arma::vec tmp_w(tmp);
  w=tmp_w;
}

void init_full(mat_t* matfp) {
  init_extended(matfp);

  arma::mat tmp=loadMatrix(matfp, (char*) "m");
  m=tmp(0, 0);

  c=cumsum(s);
}

// ---------------------------------------------------------------------------------
// This transformation function is used to break the symmetry of symmetric
// functions.
// ---------------------------------------------------------------------------------
arma::mat T_asy(const arma::mat& f, long double beta) {
  arma::mat g=f;
  arma::mat temp=beta*(arma::linspace(0, 1, f.n_rows).t());
  arma::uvec ind=find(f > 0);

  for (unsigned i=0; i<ind.n_rows; i++)
    g(ind(i))=pow(f(ind(i)), 1+(temp(ind(i)) * sqrt(f(ind(i)))));

  return g;
}

// ---------------------------------------------------------------------------------
// This transformation is used to create the ill-conditioning effect.
// ---------------------------------------------------------------------------------
arma::mat T_diag(const arma::mat& f, long double alpha) {
  long double alphasq=sqrt(alpha);

  arma::mat lin=arma::linspace(0, 1, f.n_rows);
  arma::mat linpow(1, f.n_rows);

  for (unsigned i=0; i<linpow.n_cols; i++)
    linpow(0, i)=pow(alphasq, lin(i, 0));

  // ToDo: Check if we can do this more efficient by transposing in the loop
  arma::mat scales=linpow.t();

  return scales % f;
}

// ---------------------------------------------------------------------------------
// This transformation is used to create smooth local irregularities.
//----------------------------------------------------------------------------------
arma::mat T_irreg(const arma::mat& f) {
  long double a=0.1;
  arma::mat g=f;
  arma::uvec idx=find(f > 0);
  g.elem(idx)=log(f.elem(idx))/a;
  g.elem(idx)=pow(exp(g.elem(idx)+0.49*(sin(g.elem(idx))+sin(0.79*g.elem(idx)))), a);
  idx=find(f < 0);
  g.elem(idx)=log(-f.elem(idx))/a;
  // ToDo: asegurarme de quien tiene mas prioridad, si pow o el signo menos
  g.elem(idx)=-pow(exp(g.elem(idx)+0.49*(sin(0.55*g.elem(idx))+sin(0.31*g.elem(idx)))), a);
  return g;
}

// ---------------------------------------------------------------------------------
// This function tests a given decision vector against the boundaries of a function.
// ---------------------------------------------------------------------------------
arma::uvec checkBounds(const arma::mat& x, long double lb, long double ub) {
  arma::uvec idx1=find(x > ub);
  arma::uvec idx2=find(x < lb);
  return join_cols(idx1, idx2);
}

// /////////////// //
// Basic functions //
// /////////////// //

// ---------------------------------------------------------------------------------
// Sphere Function
// ---------------------------------------------------------------------------------
long double sphere(const arma::mat& x) {
  arma::mat fit=sum(square(x), 0);
  return fit(0);
}

// ---------------------------------------------------------------------------------
// Elliptic Function
// ---------------------------------------------------------------------------------
long double elliptic(const arma::mat& x) {
  long double condition=1e+06;

  arma::mat lin=arma::linspace(0, 1, x.n_rows);
  arma::mat coefficients(1, x.n_rows);

  for (unsigned i=0; i<coefficients.n_cols; i++)
    coefficients(0, i)=pow(condition, lin(i, 0));

  long double fit=dot(coefficients, square(T_irreg(x)));
  return fit;
}

// ---------------------------------------------------------------------------------
// Rastrigin's Function
// ---------------------------------------------------------------------------------
long double rastrigin(const arma::mat& x) {
  long double A=10;
  arma::mat newX=T_diag(T_asy(T_irreg(x), 0.2), 10);
  arma::mat fit=A*(x.n_rows-sum(cos(2*arma::datum::pi*newX), 0)) + sum(square(newX), 0);
  return fit(0);
}

// ---------------------------------------------------------------------------------
// Ackley's Function
// ---------------------------------------------------------------------------------
long double ackley(const arma::mat& x) {
  arma::mat tmpsq =sum(square(x), 0);
  arma::mat tmpcos=sum(cos(2*arma::datum::pi*x), 0);
  long double fit=20-20*exp(-0.2*sqrt(tmpsq(0)/x.n_rows))-exp(tmpcos(0)/x.n_rows)+exp(1);
  return fit;
}

// ---------------------------------------------------------------------------------
// Schwefel's Problem 1.2
// ---------------------------------------------------------------------------------
long double schwefel(const arma::mat& x) {
  arma::mat newX=T_asy(T_irreg(x), 0.2);
  long double fit=0;

  for (unsigned i=0; i<newX.n_rows; i++) {
    arma::mat f=sum(newX.submat(0, 0, i, 0));
    fit=fit+pow(f(0), 2);
  }

  return fit;
}

// ---------------------------------------------------------------------------------
// Rosenbrock's Function
// ---------------------------------------------------------------------------------
long double rosenbrock(const arma::mat& x) {
  arma::mat fit=sum(100*(square(square(x.submat(0, 0, x.n_rows-2, 0))-x.submat(1, 0, x.n_rows-1, 0)))+square(x.submat(0, 0, x.n_rows-2, 0)-1));
  return fit(0);
}

// /////////////////// //
// Benchmark functions //
// /////////////////// //

// ---------------------------------------------------------------------------------
// f1: Shifted Elliptic Function
// D = 1000
// ---------------------------------------------------------------------------------
void init_f1(mat_t* matfp) {
  init_simple(matfp);
}

long double f1(const arma::mat& x) {
  arma::uvec idx=checkBounds(x, lb, ub);

  // We do not allow evaluating more than one solution at the same time
  arma::mat newX=x-xopt;

  long double fit=elliptic(newX);

  if (idx.n_rows*idx.n_cols > 0) {
    fit=std::numeric_limits<long double>::max();
    std::cerr << "Some of the solutions are violating boundary constraints." << std::endl;
  }

  return fit;
}

// ---------------------------------------------------------------------------------
// f2: Shifted Rastrigin's Function
// D = 1000
// ---------------------------------------------------------------------------------
void init_f2(mat_t* matfp) {
  init_simple(matfp);
}

long double f2(const arma::mat& x) {
  arma::uvec idx=checkBounds(x, lb, ub);

  // We do not allow evaluating more than one solution at the same time
  arma::mat newX=x-xopt;

  long double fit=rastrigin(newX);

  if (idx.n_rows*idx.n_cols > 0) {
    fit=std::numeric_limits<long double>::max();
    std::cerr << "Some of the solutions are violating boundary constraints." << std::endl;
  }

  return fit;
}

// ---------------------------------------------------------------------------------
// f3: Shifted Ackley's Function
// D = 1000
// ---------------------------------------------------------------------------------
void init_f3(mat_t* matfp) {
  init_simple(matfp);
}

long double f3(const arma::mat& x) {
  arma::uvec idx=checkBounds(x, lb, ub);

  // We do not allow evaluating more than one solution at the same time
  arma::mat newX=x-xopt;

  long double fit=ackley(newX);

  if (idx.n_rows*idx.n_cols > 0) {
    fit=std::numeric_limits<long double>::max();
    std::cerr << "Some of the solutions are violating boundary constraints." << std::endl;
  }

  return fit;
}

// ---------------------------------------------------------------------------------
// f4: 7-nonseparable, 1-separable Shifted and Rotated Elliptic Function
// D = 1000
// ---------------------------------------------------------------------------------
void init_f4(mat_t* matfp) {
  init_extended(matfp);
}

long double f4(const arma::mat& x) {
  arma::uvec idx=checkBounds(x, lb, ub);

  // We do not allow evaluating more than one solution at the same time
  arma::mat newX=x-xopt;

  long double fit=0.0;
  unsigned ldim=0;

  for (unsigned i=0; i<s.n_rows; i++) {
    long double f=0;
    if      (s(i)==25 ) {
      f=elliptic(R25*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    else if (s(i)==50 ) {
      f=elliptic(R50*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    else if (s(i)==100) {
      f=elliptic(R100*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    fit=fit+w(i)*f;
  }

  fit=fit+elliptic(newX.elem(p.subvec(ldim, p.n_rows-1)));

  if (idx.n_rows*idx.n_cols > 0) {
    fit=std::numeric_limits<long double>::max();
    std::cerr << "Some of the solutions are violating boundary constraints." << std::endl;
  }

  return fit;
}

// ---------------------------------------------------------------------------------
// f5: 7-nonseparable, 1-separable Shifted and Rotated Rastrigin's Function
// D = 1000
// ---------------------------------------------------------------------------------
void init_f5(mat_t* matfp) {
  init_extended(matfp);
}

long double f5(const arma::mat& x) {
  arma::uvec idx=checkBounds(x, lb, ub);

  // We do not allow evaluating more than one solution at the same time
  arma::mat newX=x-xopt;

  long double fit=0.0;
  unsigned ldim=0;

  for (unsigned i=0; i<s.n_rows; i++) {
    long double f=0;
    if      (s(i)==25 ) {
      f=rastrigin(R25*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    else if (s(i)==50 ) {
      f=rastrigin(R50*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    else if (s(i)==100) {
      f=rastrigin(R100*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    fit=fit+w(i)*f;
  }

  fit=fit+rastrigin(newX.elem(p.subvec(ldim, p.n_rows-1)));

  if (idx.n_rows*idx.n_cols > 0) {
    fit=std::numeric_limits<long double>::max();
    std::cerr << "Some of the solutions are violating boundary constraints." << std::endl;
  }

  return fit;
}

// ---------------------------------------------------------------------------------
// f6: 7-nonseparable, 1-separable Shifted and Rotated Ackley's Function
// D = 1000
// ---------------------------------------------------------------------------------
void init_f6(mat_t* matfp) {
  init_extended(matfp);
}

long double f6(const arma::mat& x) {
  arma::uvec idx=checkBounds(x, lb, ub);

  // We do not allow evaluating more than one solution at the same time
  arma::mat newX=x-xopt;

  long double fit=0.0;
  unsigned ldim=0;

  for (unsigned i=0; i<s.n_rows; i++) {
    long double f=0;
    if      (s(i)==25 ) {
      f=ackley(R25*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    else if (s(i)==50 ) {
      f=ackley(R50*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    else if (s(i)==100) {
      f=ackley(R100*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    fit=fit+w(i)*f;
  }

  fit=fit+ackley(newX.elem(p.subvec(ldim, p.n_rows-1)));

  if (idx.n_rows*idx.n_cols > 0) {
    fit=std::numeric_limits<long double>::max();
    std::cerr << "Some of the solutions are violating boundary constraints." << std::endl;
  }

  return fit;
}

// ---------------------------------------------------------------------------------
// f7: 7-nonseparable, 1-separable Shifted Schwefel's Function
// D = 1000
// ---------------------------------------------------------------------------------
void init_f7(mat_t* matfp) {
  init_extended(matfp);
}

long double f7(const arma::mat& x) {
  arma::uvec idx=checkBounds(x, lb, ub);

  // We do not allow evaluating more than one solution at the same time
  arma::mat newX=x-xopt;

  long double fit=0.0;
  unsigned ldim=0;

  for (unsigned i=0; i<s.n_rows; i++) {
    long double f=0;
    if      (s(i)==25 ) {
      f=schwefel(R25*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    else if (s(i)==50 ) {
      f=schwefel(R50*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    else if (s(i)==100) {
      f=schwefel(R100*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    fit=fit+w(i)*f;
  }

  fit=fit+sphere(newX.elem(p.subvec(ldim, p.n_rows-1)));

  if (idx.n_rows*idx.n_cols > 0) {
    fit=std::numeric_limits<long double>::max();
    std::cerr << "Some of the solutions are violating boundary constraints." << std::endl;
  }

  return fit;
}

// ---------------------------------------------------------------------------------
// f8: 20-nonseparable Shifted and Rotated Elliptic Function
// D = 1000
// ---------------------------------------------------------------------------------
void init_f8(mat_t* matfp) {
  init_extended(matfp);
}

long double f8(const arma::mat& x) {
  arma::uvec idx=checkBounds(x, lb, ub);

  // We do not allow evaluating more than one solution at the same time
  arma::mat newX=x-xopt;

  long double fit=0.0;
  unsigned ldim=0;

  for (unsigned i=0; i<s.n_rows; i++) {
    long double f=0;
    if      (s(i)==25 ) {
      f=elliptic(R25*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    else if (s(i)==50 ) {
      f=elliptic(R50*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    else if (s(i)==100) {
      f=elliptic(R100*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    fit=fit+w(i)*f;
  }

  if (idx.n_rows*idx.n_cols > 0) {
    fit=std::numeric_limits<long double>::max();
    std::cerr << "Some of the solutions are violating boundary constraints." << std::endl;
  }

  return fit;
}

// ---------------------------------------------------------------------------------
// f9: 20-nonseparable Shifted and Rotated Rastrigin's Function
// D = 1000
// ---------------------------------------------------------------------------------
void init_f9(mat_t* matfp) {
  init_extended(matfp);
}

long double f9(const arma::mat& x) {
  arma::uvec idx=checkBounds(x, lb, ub);

  // We do not allow evaluating more than one solution at the same time
  arma::mat newX=x-xopt;

  long double fit=0.0;
  unsigned ldim=0;

  for (unsigned i=0; i<s.n_rows; i++) {
    long double f=0;
    if      (s(i)==25 ) {
      f=rastrigin(R25*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    else if (s(i)==50 ) {
      f=rastrigin(R50*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    else if (s(i)==100) {
      f=rastrigin(R100*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    fit=fit+w(i)*f;
  }

  if (idx.n_rows*idx.n_cols > 0) {
    fit=std::numeric_limits<long double>::max();
    std::cerr << "Some of the solutions are violating boundary constraints." << std::endl;
  }

  return fit;
}

// ---------------------------------------------------------------------------------
// f10: 20-nonseparable Shifted and Rotated Ackley's Function
// D = 1000
// ---------------------------------------------------------------------------------
void init_f10(mat_t* matfp) {
  init_extended(matfp);
}

long double f10(const arma::mat& x) {
  arma::uvec idx=checkBounds(x, lb, ub);

  // We do not allow evaluating more than one solution at the same time
  arma::mat newX=x-xopt;

  long double fit=0.0;
  unsigned ldim=0;

  for (unsigned i=0; i<s.n_rows; i++) {
    long double f=0;
    if      (s(i)==25 ) {
      f=ackley(R25*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    else if (s(i)==50 ) {
      f=ackley(R50*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    else if (s(i)==100) {
      f=ackley(R100*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    fit=fit+w(i)*f;
  }

  if (idx.n_rows*idx.n_cols > 0) {
    fit=std::numeric_limits<long double>::max();
    std::cerr << "Some of the solutions are violating boundary constraints." << std::endl;
  }

  return fit;
}

// ---------------------------------------------------------------------------------
// f11: 20-nonseparable Shifted Schwefel's Function
// D = 1000
// ---------------------------------------------------------------------------------
void init_f11(mat_t* matfp) {
  init_extended(matfp);
}

long double f11(const arma::mat& x) {
  arma::uvec idx=checkBounds(x, lb, ub);

  // We do not allow evaluating more than one solution at the same time
  arma::mat newX=x-xopt;

  long double fit=0.0;
  unsigned ldim=0;

  for (unsigned i=0; i<s.n_rows; i++) {
    long double f=0;
    if      (s(i)==25 ) {
      f=schwefel(R25*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    else if (s(i)==50 ) {
      f=schwefel(R50*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    else if (s(i)==100) {
      f=schwefel(R100*newX.elem(p.subvec(ldim, ldim+s(i)-1)));
      ldim = ldim + s(i);
    }
    fit=fit+w(i)*f;
  }

  if (idx.n_rows*idx.n_cols > 0) {
    fit=std::numeric_limits<long double>::max();
    std::cerr << "Some of the solutions are violating boundary constraints." << std::endl;
  }

  return fit;
}

// ---------------------------------------------------------------------------------
// f12: Shifted Rosenbrock's Function
// D = 1000
// ---------------------------------------------------------------------------------
void init_f12(mat_t* matfp) {
  init_simple(matfp);
}

long double f12(const arma::mat& x) {
  arma::uvec idx=checkBounds(x, lb, ub);

  // We do not allow evaluating more than one solution at the same time
  arma::mat newX=x-xopt;

  long double fit=rosenbrock(newX);

  if (idx.n_rows*idx.n_cols > 0) {
    fit=std::numeric_limits<long double>::max();
    std::cerr << "Some of the solutions are violating boundary constraints." << std::endl;
  }

  return fit;
}

// ---------------------------------------------------------------------------------
// f13: Shifted Schwefel's Function with Conforming Overlapping Subcomponents
// D = 905
// ---------------------------------------------------------------------------------
void init_f13(mat_t* matfp) {
  init_full(matfp);
}

long double f13(const arma::mat& x) {
  arma::uvec idx=checkBounds(x, lb, ub);

  // We do not allow evaluating more than one solution at the same time
  arma::mat newX=x-xopt;

  long double fit=0.0;

  for (unsigned i=0; i<s.n_rows; i++) {
    unsigned ldim=0;
    long double f=0;

    if (i==0)
      ldim=0;
    else
      ldim=c(i-1)-(i*m);

    unsigned udim=c(i)-(i*m)-1;

    if      (s(i)==25 )
      f=schwefel(R25*newX.elem(p.subvec(ldim, udim)));
    else if (s(i)==50 )
      f=schwefel(R50*newX.elem(p.subvec(ldim, udim)));
    else if (s(i)==100)
      f=schwefel(R100*newX.elem(p.subvec(ldim, udim)));

    fit=fit+w(i)*f;
  }

  if (idx.n_rows*idx.n_cols > 0) {
    fit=std::numeric_limits<long double>::max();
    std::cerr << "Some of the solutions are violating boundary constraints." << std::endl;
  }

  return fit;
}

// ---------------------------------------------------------------------------------
// f14: Shifted Schwefel's Function with Conflicting Overlapping Subcomponents
// D = 905
// ---------------------------------------------------------------------------------
void init_f14(mat_t* matfp) {
  init_full(matfp);
}

long double f14(const arma::mat& x) {
  arma::uvec idx=checkBounds(x, lb, ub);

  long double fit=0.0;

  for (unsigned i=0; i<s.n_rows; i++) {
    unsigned ldim=0, ldimshift=0;
    long double f=0;

    if (i==0) {
      ldim=0;
      ldimshift=0;
    }
    else {
      ldim=c(i-1)-(i*m);
      ldimshift=c(i-1);
    }

    unsigned udim=c(i)-(i*m)-1;
    unsigned udimshift=c(i)-1;

    arma::mat z=(x.elem(p.subvec(ldim, udim))-xopt.submat(ldimshift, 0, udimshift, 0));

    if      (s(i)==25 )
      f=schwefel(R25*z);
    else if (s(i)==50 )
      f=schwefel(R50*z);
    else if (s(i)==100)
      f=schwefel(R100*z);

    fit=fit+w(i)*f;
  }

  if (idx.n_rows*idx.n_cols > 0) {
    fit=std::numeric_limits<long double>::max();
    std::cerr << "Some of the solutions are violating boundary constraints." << std::endl;
  }

  return fit;
}

// ---------------------------------------------------------------------------------
// f15: Shifted Schwefel's Function
// D = 1000
// ---------------------------------------------------------------------------------
void init_f15(mat_t* matfp) {
  init_simple(matfp);
}

long double f15(const arma::mat& x) {
  arma::uvec idx=checkBounds(x, lb, ub);

  // We do not allow evaluating more than one solution at the same time
  arma::mat newX=x-xopt;

  long double fit=schwefel(newX);

  if (idx.n_rows*idx.n_cols > 0) {
    fit=std::numeric_limits<long double>::max();
    std::cerr << "Some of the solutions are violating boundary constraints." << std::endl;
  }

  return fit;
}
