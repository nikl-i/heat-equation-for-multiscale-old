/*
 *	Tridiagonal matrix algorithm.
 *	307 group, Nikishin Nikolay.
 * 10^(-15)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>

#define EPSILON 0.00000001

double *
tdma (const int n, const  double *C, const double *A, const double *B, const double *F,
		const double k1, const double m1, const double k2, const double m2)
{
	int i;
	double *x,*alpha,*beta,temp;

	alpha = (double *) calloc(n,sizeof(double));
	if(alpha == NULL)
	{
		fprintf(stderr,"Error: Memory allocation Error.\n");
		perror("calloc");
		return NULL;
	}
	beta = (double *) calloc(n,sizeof(double));
	if(beta == NULL)
	{
		fprintf(stderr,"Error: Memory allocation Error.\n");
		perror("calloc");
		free(alpha);
		return NULL;
	}
	alpha--; /*indices starts from 1*/
	beta--;
	alpha[1]=k1;
	beta[1]=m1;

	for(i=1;i<=n-1;i++)
	{/*The first pass (setting coefficients)*/
		if(	!(abs(C[i]) >= abs(A[i]) + abs(B[i])))
			fprintf(stderr,"Warning! Violation of diagonal supremacy condition at %d string.\n", i);
		temp = C[i] - A[i]*alpha[i];
		/*printf("%d| A[i] %lf C[i] %lf alpha[i] %lf temp %lf \n",i,A[i],C[i],alpha[i],temp);*/
		if(fabs(temp) < EPSILON)
		{
			free(++alpha);
			free(++beta);
			fprintf(stderr,"Error: Division by zero at calculation %d coefficients.\n", i);
			return NULL;
		}
		/*printf("%d| B[i] %lf A[i] %lf F[i] %lf beta[i] %lf temp %lf\n",i,B[i],A[i],F[i],beta[i],temp);*/
		alpha[i+1] = B[i] / temp;
		beta[i+1] = ( F[i] + A[i] * beta[i] ) / temp;
		/*printf("%d| a %lf  b %lf\n",i,alpha[i+1],beta[i+1]);*/
	}

	x = (double *) calloc(n+1,sizeof(double));

	if(x == NULL)
	{
		fprintf(stderr,"Error: Memory allocation Error.\n");
		perror("calloc");
		free(++alpha);
		free(++beta);
		return NULL;
	}
	x[n] = (m2 + k2 * beta[n] ) / ( 1 - alpha[n] * k2);
	/*printf("m2 %lf k2 %lf beta_n %lf alpha_n %lf | x: %lf \n",m2,k2,beta[n],alpha[n],x[n]);*/
	
	for (i=n; i>=1; i--)
	{/* The second pass (back-substition) */
		x[i-1]=alpha[i]*x[i] + beta[i];
	}

	free(++alpha);
	free(++beta);
	return x;
}

void
mem_free(double *x,double *C,double *A,double *B,double *F)
{
	if(C != NULL)
		free(C);
	if(B != NULL)
		free(B);
	if(A != NULL)
		free(A);
	if(F != NULL)
		free(F);
	if(x != NULL)
		free(x);
	return;
}

int
main (void)
{
	int i,n,check;
	double *x=NULL,*C=NULL,*A=NULL,*B=NULL,*F=NULL;
	double p0,q0,pn,qn,max,temp=0;

	const double pi=4*atan(1);

/*	
	double test[]={0.059762610893402,-0.021619301791879,-0.41661856640289,-1.021737073764143,8.300883066932958};
*/

/*	
	double test[]={4.166666666666666,3.166666666666666,4.333333333333333,5.75,6.75};
*/

/*	fprintf(stdout,"Input n:\n");*/
	check=scanf("%d",&n);
	if(check<=0)
	{
		fprintf(stderr,"Error: Input Error.\n");
		perror("scanf");
		return EXIT_FAILURE;
	}
	if (!n)
		return 0;

	B = (double *) calloc(n+1,sizeof(double));
	if(B == NULL)
	{
		fprintf(stderr,"Error: Memory allocation Error.\n");
		perror("calloc");
		mem_free(x,C,A,B,F);
		return EXIT_FAILURE;
	}

/*	fprintf(stdout,"Input diagonal elements:\n");*/
	for (i=0;i<=n-2;i++)
	{
		check=scanf("%lf",&B[i]);
		if(check<=0)
		{
			fprintf(stderr,"Error: Input Error.\n");
			perror("scanf");
			mem_free(x,C,A,B,F);
			return EXIT_FAILURE;
		}
		/*printf("B[%d]: %lf\n",i,B[i]);*/
	}
	
	A = (double *) calloc(n+1,sizeof(double));
	if(A == NULL)
	{
		mem_free(x,C,A,B,F);
		fprintf(stderr,"Error: Memory allocation Error.\n");
		perror("calloc");
		return EXIT_FAILURE;
	}

/*	fprintf(stdout,"Input below diagonal elements:\n");*/
	for (i=0;i<=n-2;i++)
	{
		check=scanf("%lf",&A[i]);
		if(check<=0)
		{
			fprintf(stderr,"Error: Input Error.\n");
			perror("scanf");
			mem_free(x,C,A,B,F);
			return EXIT_FAILURE;
		}
	}

	C = (double *) calloc(n+1,sizeof(double));
	if(C == NULL)
	{
		fprintf(stderr,"Error: Memory allocation Error.\n");
		perror("calloc");
		return EXIT_FAILURE;
	}

/*	fprintf(stdout,"Input diagonal elements:\n");*/
	for (i=0;i<=n-2;i++)
	{
		check=scanf("%lf",&C[i]);
		if(check<=0)
		{
			fprintf(stderr,"Error: Input Error.\n");
			perror("scanf");
			mem_free(x,C,A,B,F);
			return EXIT_FAILURE;
		}
	}
	
	F = (double *) calloc(n+1,sizeof(double));
	if(F == NULL)
	{
		fprintf(stderr,"Error: Memory allocation Error.\n");
		perror("calloc");
		mem_free(x,C,A,B,F);
		return EXIT_FAILURE;
	}

/*	fprintf(stdout,"Введите элементы правой части:\n");*/
	for (i=0;i<=n-2;i++)
	{
		check=scanf("%lf",&F[i]);
		if(check<=0)
		{
			mem_free(x,C,A,B,F);
			fprintf(stderr,"Error: Input Error.\n");
			perror("scanf");
			return EXIT_FAILURE;
		}
	}

/*	fprintf(stdout,"Введите левое граничное условие:\n");*/
	check=scanf("%lf %lf",&p0,&q0);
	if(check<=0)
	{
		mem_free(x,C,A,B,F);
		fprintf(stderr,"Error: Input Error.\n");
		perror("scanf");
		return EXIT_FAILURE;
	}

/*	fprintf(stdout,"Введите правое граничное условие:\n");*/
	check=scanf("%lf %lf",&pn,&qn);
	if(check<=0)
	{
		mem_free(x,C,A,B,F);
		fprintf(stderr,"Error: Input Error.\n");
		perror("scanf");
		return EXIT_FAILURE;
	}
	
	/* A, B, C, F indices starts from 1 !*/
	x=tdma(n,C-1,A-1,B-1,F-1,p0,q0,pn,qn);

	for(i=0;i<n+1;i++)
		printf("%.64lf\n",x[i]);
	printf("\n");
	
	max=fabs(x[0]);
	for(i=1;i<=n;i++)
	{
		temp=fabs(sin(pi*i/n) - x[i]);
		printf("x[%d]: %.64lf\ny[%d]: %.64lf\nmisalignment[%d]: %.64lf\n",i,x[i],i,(double)sin(pi*i/n),i,temp);
		if (max<temp) 
			max=temp;
	}
		
	/*
	for(i=0;i<=n1;i++)
	{
		temp=fabs(test[i] - x[i]);
		printf("%d\n x: %.64lf\n y: %.64lf\n diff: %.64lf \n ",i,x[i],test[i], temp);
		if (max<temp) 
			max=temp;
	}
	*/
	printf("Calculating Error: %.64lf\n",max);
	mem_free(x,C,A,B,F);
	
	return EXIT_SUCCESS;
}
