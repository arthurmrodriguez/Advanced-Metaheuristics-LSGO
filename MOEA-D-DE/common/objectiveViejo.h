#ifndef __OBJECTIVE_H_
#define __OBJECTIVE_H_

#include "../common/global.h"

// control the PF shape
void alphafunction(double alpha[], vector<double>&x, int dim, int type)
{
    if(dim==2)
	{
        if(type==21){
		    alpha[0] = x[0];
			alpha[1] = 1 - sqrt(x[0]);
		}

        if(type==22){
		    alpha[0] = x[0];
			alpha[1] = 1 - x[0]*x[0];
		}

        if(type==23){
		    alpha[0] = x[0];
			alpha[1] = 1 - sqrt(alpha[0]) - alpha[0]*sin(10*alpha[0]*alpha[0]*pi);
		}

		if(type==24){
		    alpha[0] = x[0];
			alpha[1] = 1 - x[0] - 0.05*sin(4*pi*x[0]);
		}

	}
	else
	{

		if(type==31){
			alpha[0] = cos(x[0]*pi/2)*cos(x[1]*pi/2);
			alpha[1] = cos(x[0]*pi/2)*sin(x[1]*pi/2);
			alpha[2] = sin(x[0]*pi/2);
		}

		if(type==32){
			alpha[0] = 1 - cos(x[0]*pi/2)*cos(x[1]*pi/2);
			alpha[1] = 1 - cos(x[0]*pi/2)*sin(x[1]*pi/2);
			alpha[2] = 1 - sin(x[0]*pi/2);
		}

		if(type==33){
		    alpha[0] = x[0];
			alpha[1] = x[1];
			alpha[2] = 3 - (sin(3*pi*x[0]) + sin(3*pi*x[1])) - 2*(x[0] + x[1]);
		}

		if(type==34){
			alpha[0] = x[0]*x[1];
			alpha[1] = x[0]*(1 - x[1]);
			alpha[2] = (1 - x[0]);
		}
	}
}

// control the distance
double betafunction(vector <double>x, int type)
{
	double beta;
	int dim = x.size();

	if(dim==0) beta = 0;

    if(type==1){
		beta = 0;
		for(int i=0; i<dim; i++){
		    beta+= x[i]*x[i];
		}
		beta = 2.0*beta/dim;
	}

    if(type==2){
		beta = 0;
		for(int i=0; i<dim; i++){
		    beta+= sqrt(i+1)*x[i]*x[i];
		}
		beta = 2.0*beta/dim;
	}

	if(type==3){
		double sum = 0, xx;
		for(int i=0; i<dim; i++){
			xx = 2*x[i];
		    sum+= (xx*xx - cos(4*pi*xx) + 1);
		}
	    beta = 2.0*sum/dim;
	}

	if(type==4){
		double sum = 0, prod = 1, xx;
		for(int i=0; i<dim; i++){
			xx  = 2*x[i];
		    sum+= xx*xx;
			prod*=cos(10*pi*xx/sqrt(i+1));
		}
		beta = 2.0*(sum - 2*prod + 2)/dim;
	}

	return beta;
}


// control the PS shape of 2-d instances
double psfunc2(double &x, double &t1, int dim, int type, int css){
	// type:  the type of curve
	// css:   the class of index
	double beta;

	dim++;

	if(type==21){
		double xy   = 2*(x - 0.5);
		beta = xy - pow(t1, 0.5*(nvar + 3*dim - 8)/(nvar - 2));
	}

	if(type==22){
		double theta = 6*pi*t1 + dim*pi/nvar;
		double xy    = 2*(x - 0.5);
		beta = xy - sin(theta);
	}

	if(type==23){
		double theta = 6*pi*t1 + dim*pi/nvar;
		double ra    = 0.8*t1;
		double xy    = 2*(x - 0.5);
		if(css==1)
			beta = xy - ra*cos(theta);
		else{
			beta = xy - ra*sin(theta);
		}
	}

	if(type==24){
		double theta = 6*pi*t1 + dim*pi/nvar;
		double xy    = 2*(x - 0.5);
		double ra    = 0.8*t1;
		if(css==1)
			beta = xy - ra*cos(theta/3);
		else{
			beta = xy - ra*sin(theta);
		}
	}

	if(type==25){
        double rho   = 0.8;
		double phi   = pi*t1;
		double theta = 6*pi*t1 + dim*pi/nvar;
		double xy    = 2*(x - 0.5);
		if(css==1)
			beta = xy - rho*sin(phi)*sin(theta);
		else if(css==2)
			beta = xy - rho*sin(phi)*cos(theta);
		else
			beta = xy - rho*cos(phi);
	}

	if(type==26){
		double theta = 6*pi*t1 + dim*pi/nvar;
		double ra    = 0.3*t1*(t1*cos(4*theta) + 2);
		double xy    = 2*(x - 0.5);
		if(css==1)
			beta = xy - ra*cos(theta);
		else{
			beta = xy - ra*sin(theta);
		}
	}

	return beta;
}


// control the PS shapes of 3-D instances
double psfunc3(double &x, double &t1, double &t2, int dim, int type){
	// type:  the type of curve
	// css:   the class of index
	double beta;

	dim++;

	if(type==31){
		double xy  = 4*(x - 0.5);
		double rate = 1.0*dim/nvar;
		beta = xy - 4*(t1*t1*rate + t2*(1.0-rate)) + 2;
	}

	if(type==32){
		double theta = 2*pi*t1 + dim*pi/nvar;
		double xy    = 4*(x - 0.5);
		beta = xy - 2*t2*sin(theta);
	}

	return beta;
}

void objective(vector<double> &x_var, vector <double> &y_obj)
{
	// 2-objective case
	if(nobj==2)
	{
		if(ltype==21||ltype==22||ltype==23||ltype==24||ltype==26)
		{
			double g = 0, h = 0, a, b;
			vector <double> aa;
			vector <double> bb;
			for(int n=1;n<nvar;n++)
			{

				if(n%2==0){
					a = psfunc2(x_var[n],x_var[0],n,ltype,1);  // linkage
					aa.push_back(a);
				}
				else
				{
					b = psfunc2(x_var[n],x_var[0],n,ltype,2);
					bb.push_back(b);
				}

			}

			g = betafunction(aa,dtype);
			h = betafunction(bb,dtype);

			double alpha[2];
			alphafunction(alpha,x_var,2,ptype);  // shape function
			y_obj[0] = alpha[0] + h;
			y_obj[1] = alpha[1] + g;
			aa.clear();
			bb.clear();
		}

		if(ltype==25)
		{
			double g = 0, h = 0, a, b;
			double e = 0, c;
			vector <double> aa;
			vector <double> bb;
			for(int n=1;n<nvar;n++){
				if(n%3==0){
					a = psfunc2(x_var[n],x_var[0],n,ltype,1);
					aa.push_back(a);
				}
				else if(n%3==1)
				{
					b = psfunc2(x_var[n],x_var[0],n,ltype,2);
					bb.push_back(b);
				}
				else{
					c = psfunc2(x_var[n],x_var[0],n,ltype,3);
					if(n%2==0)    aa.push_back(c);
					else          bb.push_back(c);
				}
			}
			g = betafunction(aa,dtype);          // distance function
			h = betafunction(bb,dtype);
			double alpha[2];
			alphafunction(alpha,x_var,2,ptype);  // shape function
			y_obj[0] = alpha[0] + h;
			y_obj[1] = alpha[1] + g;
			aa.clear();
			bb.clear();
		}
	}


	// 3-objective case
	if(nobj==3)
	{
		if(ltype==31||ltype==32)
		{
			double g = 0, h = 0, e = 0, a;
			vector <double> aa;
			vector <double> bb;
			vector <double> cc;
			for(int n=2;n<nvar;n++)
			{
				a = psfunc3(x_var[n],x_var[0],x_var[1],n,ltype);
				if(n%3==0)	    aa.push_back(a);
				else if(n%3==1)	bb.push_back(a);
				else            cc.push_back(a);
			}

			g = betafunction(aa,dtype);
			h = betafunction(bb,dtype);
			e = betafunction(cc,dtype);

			double alpha[3];
			alphafunction(alpha,x_var,3,ptype);  // shape function
			y_obj[0] = alpha[0] + h;
			y_obj[1] = alpha[1] + g;
			y_obj[2] = alpha[2] + e;
			aa.clear();
			bb.clear();
			cc.clear();
		}
	}
}


#endif
