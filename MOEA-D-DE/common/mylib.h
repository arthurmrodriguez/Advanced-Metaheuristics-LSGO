#ifndef __MYLIB_H_
#define __MYLIB_H_

void minfastsort(double x[], int idx[], int n, int m)
{
    for(int i=0; i<m; i++)
	{
	    for(int j=i+1; j<n; j++)
			if(x[i]>x[j])
			{
			    double temp = x[i];
				x[i]        = x[j];
				x[j]        = temp;
				int id      = idx[i];
				idx[i]      = idx[j];
				idx[j]      = id;
			}
	}
}

double dist_array(double vec1[], double vec2[], int dim)
{
    double sum = 0;
	for(int n=0; n<dim; n++)
	    sum+= (vec1[n] - vec2[n])*(vec1[n] - vec2[n]);
	return sqrt(sum);
}

double dist_vector(vector <double> &vec1, vector <double> &vec2)
{
	int dim = vec1.size();
    double sum = 0;
	for(int n=0; n<dim; n++)
	    sum+=(vec1[n] - vec2[n])*(vec1[n] - vec2[n]);
	return sqrt(sum);
}

double norm_vector(vector <double> &vec)
{
	int    dim = vec.size();
	double sum = 0;
	for(int i=0;i<dim;i++)
        sum = sum + vec[i]*vec[i];
    return sqrt(sum);
}



double sum_vector(vector<double>&vec)
{
	int dim = vec.size();
	double sum = 0;
	for(int i=0;i<dim;i++)
        sum = sum + vec[i];
    return sum;
}

// inner product
double prod_vector(vector <double>&vec1, vector <double>&vec2)
{
	int dim = vec1.size();
    double sum = 0;
	for(int i=0; i<dim; i++)
		sum+= vec1[i]*vec2[i];
	return sum;
}

#endif
