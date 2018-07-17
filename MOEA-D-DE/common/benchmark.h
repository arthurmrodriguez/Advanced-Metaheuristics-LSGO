#ifndef _BENCHMARK_H_
#define _BENCHMARK_H_

#define pi 3.1415926

void benchmark(vector<double> x_var, vector <double> &y_obj)
{

    if(!strcmp(strTestInstance,"F1"))
	{
		double g = 0, xy;
		for(int n=1;n<nvar;n++){
            xy  = x_var[0] - x_var[n];
			g  += xy*xy;
		}
		g = 1 + 9*g/(nvar-1);
		y_obj[0] = x_var[0];
		y_obj[1] = g*(1 - sqrt(y_obj[0]/g));
	}

	if(!strcmp(strTestInstance,"F2"))
	{
		double g = 0, xy;
		for(int n=1;n<nvar;n++){
            xy  = x_var[0] - x_var[n];
			g  += xy*xy;
		}
		g = 1 + 9*g/(nvar-1);
		y_obj[0] = sqrt(x_var[0]);
		y_obj[1] = g*(1 - pow(y_obj[0]/g,2));
	}

	if(!strcmp(strTestInstance,"F3"))
	{
		double g = 0, xy;
		for(int n=1;n<nvar;n++){
			xy  = x_var[0] - x_var[n];
			g  += xy*xy;
		}		
		g = 1 + 9*pow(g/(nvar-1),0.25);

		y_obj[0] = 1 - exp(-4*x_var[0])*pow(sin(6*pi*x_var[0]),6);
		y_obj[1] = g*(1- pow(y_obj[0]/g,2));
	}

	if(!strcmp(strTestInstance,"F4"))
	{
		double g = 0,xy;
		for(int n=2; n<nvar;n++)
		{
			xy = x_var[0] - x_var[n];
			g+= xy*xy;
		}
		y_obj[0] = (1 + g)*cos(x_var[0]*pi/2)*cos(x_var[1]*pi/2);
		y_obj[1] = (1 + g)*cos(x_var[0]*pi/2)*sin(x_var[1]*pi/2);
		y_obj[2] = (1 + g)*sin(x_var[0]*pi/2);
	}

	if(!strcmp(strTestInstance,"F5"))
	{
		double g = 0, xy;
		for(int n=1;n<nvar;n++){
            xy  = x_var[0] - x_var[n]*x_var[n];
			g  += xy*xy;
		}
		g = 1 + 9*g/(nvar-1);
		y_obj[0] = x_var[0];
		y_obj[1] = g*(1 - sqrt(y_obj[0]/g));
	}

	if(!strcmp(strTestInstance,"F6"))
	{
		double g = 0, xy;
		for(int n=1;n<nvar;n++){
            xy  = x_var[0] - x_var[n]*x_var[n];
			g  += xy*xy;
		}
		g = 1 + 9*g/(nvar-1);
		y_obj[0] = sqrt(x_var[0]);
		y_obj[1] = g*(1 - pow(y_obj[0]/g,2));
	}

	if(!strcmp(strTestInstance,"F7"))
	{
		double g = 0, xy;
		for(int n=1;n<nvar;n++){
			xy  = x_var[0] - x_var[n]*x_var[n];
			g  += xy*xy;
		}		
		g = 1 + 9*pow(g/9,0.25);

		y_obj[0] = 1 - exp(-4*x_var[0])*pow(sin(6*pi*x_var[0]),6);
		y_obj[1] = g*(1- pow(y_obj[0]/g,2));
	}

	if(!strcmp(strTestInstance,"F8"))
	{
		double g = 0,xy;
		for(int n=2; n<nvar;n++)
		{
			xy = x_var[0] - x_var[n]*x_var[n];
			g+= xy*xy;
		}
		y_obj[0] = (1 + g)*cos(x_var[0]*pi/2)*cos(x_var[1]*pi/2);
		y_obj[1] = (1 + g)*cos(x_var[0]*pi/2)*sin(x_var[1]*pi/2);
		y_obj[2] = (1 + g)*sin(x_var[0]*pi/2);
	}

	if(!strcmp(strTestInstance,"F9"))
	{
		double g = 0, h = 1;
		for(int n=1;n<nvar;n++)
		{
			double x = 10*x_var[n];
			double y = x_var[0] - x*x;
			g+= y*y;
			h*= cos(y/n);
		}
		g = g/4000 - h + 2;
		y_obj[0] = x_var[0];
		y_obj[1] = g*(1- sqrt(y_obj[0]/g));
	}

	
	if(!strcmp(strTestInstance,"F10"))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
		{
			double x = 10*x_var[n];
			double y = x_var[0] - x*x;
			g+= y*y - 10*cos(2*pi*y);
		}
		g = 1 + 10*(nvar-1) + g;
		y_obj[0] = x_var[0];
		y_obj[1] = g*(1- sqrt(y_obj[0]/g));
	}

	if(!strcmp(strTestInstance,"ZDT1"))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
			g+= x_var[n];
		g = 1 + 9*g/(nvar-1);

		y_obj[0] = x_var[0];
		y_obj[1] = g*(1 - sqrt(y_obj[0]/g));
	}


	if(!strcmp(strTestInstance,"ZDT2"))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
			g+= x_var[n];
		g = 1 + 9*g/(nvar-1);
		y_obj[0] = x_var[0];
		y_obj[1] = g*(1 - pow(y_obj[0]/g,2));
	}


	if(!strcmp(strTestInstance,"ZDT3"))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
			g+= x_var[n];
		g = 1 + 9*g/(nvar-1);

		y_obj[0] = x_var[0];
		y_obj[1] = g*(1 - sqrt(x_var[0]/g) - x_var[0]*sin(10*pi*x_var[0])/g);
	}

	if(!strcmp(strTestInstance,"ZDT4"))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
		{
			double x = 10*(x_var[n] - 0.5);
			g+= x*x - 10*cos(4*pi*x);
		}
		g = 1 + 10*(nvar-1) + g;
		y_obj[0] = x_var[0];
		y_obj[1] = g*(1- sqrt(y_obj[0]/g));
	}

	if(!strcmp(strTestInstance,"ZDT6"))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
			g+= x_var[n]/(nvar - 1);
		g = 1 + 9*pow(g,0.25) ;

		y_obj[0] = 1 - exp(-4*x_var[0])*pow(sin(6*pi*x_var[0]),6);
		y_obj[1] = g*(1- pow(y_obj[0]/g,2));
	}


	if(!strcmp(strTestInstance,"DTLZ2"))
	{
		double g = 0;
		double xx = (x_var[0] + x_var[1])/2.0;
		for(int n=2; n<nvar;n++)				
		{
			double x = 2*(x_var[n] - 0.5);
			g = g + x*x; 
		}
		y_obj[0] = (1 + g)*cos(x_var[0]*pi/2)*cos(x_var[1]*pi/2);
		y_obj[1] = (1 + g)*cos(x_var[0]*pi/2)*sin(x_var[1]*pi/2);
		y_obj[2] = (1 + g)*sin(x_var[0]*pi/2);
	}

	if(!strcmp(strTestInstance,"DTLZ1"))
	{
		double g = 0;
		for(int n=2; n<nvar;n++)
			g = g + pow(x_var[n]-0.5,2) - cos(20*pi*(x_var[n] - 0.5));
		g = 100*(nvar- 2 + g);
		y_obj[0] = (1 + g)*x_var[0]*x_var[1];
		y_obj[1] = (1 + g)*x_var[0]*(1 - x_var[1]);
		y_obj[2] = (1 + g)*(1 - x_var[0]);
	}

}

#endif
