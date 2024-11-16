#include<stdio.h>
#include<math.h>
#include<complex.h>
#define Max_Iteration 1000
#define Tolerance 1e-6

void Matrix_Multiplication(double complex A[100][100],double complex v[100],double complex result[100],int n){
	//dimension of matrix can be changed

for(int i=0;i<n;i++){
	result[i] = 0.0+0.0*I;
for(int j=0;j<n;j++){
	result[i] += A[i][j]*v[j];
}
}}

double norm_value(double complex v[100],int n){
	double sum =0.0;
	for (int i=0;i<n;i++){
		double real = creal(v[i]);
		double imag = cimag(v[i]);
		sum+= real*real + imag*imag;
	}
	return sqrt(sum);
}
void normalise(double complex v[100],int n){
	double magnitude = norm_value(v,n);
	if(magnitude !=0){
		for (int i=0;i<n;i++){
			v[i] = v[i]/magnitude;}}}

double complex rayleigh(double complex A[100][100],double complex v[100],int n){
	double complex num = 0.0 + 0.0*I;
	double complex den = 0.0 + 0.0*I;
for (int i = 0; i < n; i++) {
        double complex temp = 0.0 + 0.0 * I;
        for (int j = 0; j < n; j++) {
            temp += A[i][j]*v[j];
        }
        num += conj(v[i]) * temp;
        den += conj(v[i]) * v[i];
    }
    return num/den;
}

// Hybrid Newton-Power Method
double complex hybrid_newton(double complex A[100][100], int n) {
    double complex v[100] = {0.0 + 0.0 * I};
    double complex Av[100] = {0.0 + 0.0 * I};
    double complex lambda = 0.0 + 0.0 * I;
    double complex lambda1 = 0.0 + 0.0 * I;

    for (int i = 0; i < n; i++) {
        v[i] = 1.0 + 0.0 * I;//initialising random vector
    }
    normalise(v, n);

    for (int i=0;i <Max_Iteration;i++) {    //iteration
        Matrix_Multiplication(A,v,Av,n);
        lambda = rayleigh(A,v,n);  
        for (int j = 0;j< n;j++) {
            v[j] = Av[j];
        }
        normalise(v,n);

        if (cabs(lambda-lambda1)<Tolerance) {
            break;
        }
        lambda1= lambda;
    }

    // Newton Refinement 
    for (int i = 0; i < Max_Iteration; i++) {
        double complex shift = lambda; 

        // Solve (A - shift*I)v = v
        for (int j = 0; j < n; j++) {
            Av[i] = 0.0 + 0.0 * I;
            for (int k=0;k < n;k++) {
                if (j==k) {
                    Av[j] += (A[j][k] - shift) * v[k];
                } else {
                    Av[j]+=A[j][k]*v[k];
                }
            }}
        lambda = rayleigh(A,v,n);
        //convergence
        if (cabs(lambda-lambda1)<Tolerance) {
            break;
        }
        lambda1=lambda;
    }
    return lambda;
}
int main() {
	int n;
	printf("Enter matrix size");
	scanf("%d",&n);
	
	double complex A[100][100] = {0.0+0.0*I};
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			double real,imag;
			printf("Element %d %d(real,img) :" ,i,j);
			scanf("%lf %lf" ,&real,&imag);
			A[i][j]=real+imag*I;}}
	double complex eigen_value = hybrid_newton(A,n);
	printf("%f + %fi",creal(eigen_value),cimag(eigen_value));
	return 0;
}
