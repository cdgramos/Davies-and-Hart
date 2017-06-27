#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double autocovariance(double H, int k){
    return (0.5 * (pow(abs(k - 1), (2 * H)) - 2 * pow(abs(k), (2 * H)) +
                  pow(abs(k + 1), (2 * H))));
}

double covariance(long i, double H) {
  if (i == 0) return 1;
  else return (pow(i-1,2*H)-2*pow(i,2*H)+pow(i+1,2*H))/2;
}

double drand(){
  return (rand()+1.0)/(RAND_MAX+1.0);
}

/* normal distribution, centered on 0, std dev 1 */
double randomNormal(){
  return sqrt(-2*log(drand())) * cos(2*3.141592*drand());
}


int dft(long int length, double realSample[], double *Rk, double *Ik, long int limit){
    long int i, j;
    double arg;
    double cosarg, sinarg;
    double *tempReal = NULL;
    double *tempIm = NULL;
    double PI2 = 6.283184;
    double imaginaryPart = 1.0;

    tempReal = calloc(length, sizeof(double));
    tempIm = calloc(length, sizeof(double));

    if (tempReal == NULL || tempIm == NULL){
        printf("Error: Allocating memory\n");
        return(0);
    }

    for(i=0; i<length; i++){

        arg = (-1.0 * PI2) * i / length;

        for(j=0; j<length; j++){
            tempReal[i] += (realSample[j] * (cos(j * arg)) - imaginaryPart * (sin(j * arg)));
            tempIm[i] += (realSample[j] * (sin(j * arg)) + imaginaryPart * (cos(j * arg)));
        }
    }

    // If we need just the first X values
    for (i=0; i<limit; i++){
        Rk[i] = tempReal[i];
        Ik[i] = tempIm[i];
    }

    free(tempReal);
    free(tempIm);
    return(1);
}


int main(){
    double H = 0.05; // Desired hurst value
    int N = pow(2,10); // Number of elements
    double L = 1.0; // Lenght of realization

    double increment = pow(increment, H);
    double scale = L / N;

    double mu = 0.0; // Parameter for the normal gaussian distribution
    double sigma = 1.0; // Parameter for the normal gaussian distribution
    double rnd = 0.0; // Pseudo-random value

    double *fgn;
    double *fgn2;
    double *row_component;
    double *reverse_component;
    double *row;
    double *w;

    double *outreal;
    double *outimag;

    int i = 0;
    int j = 0;

    int flag = 1;

    FILE *FP;
    FP = fopen("file.txt", "w");

    printf("Computing a self-similar sequence with:\n");
    printf("\tHurst value: %.2lf\n\tNumber of elements: %d\n\tLength of realization: %.2lf\n",H,N,L);

    // Parameter validation (H, N, L)
    if(H <= 0 || H >= 1){
        printf("Error: H must be between ]0, 1[...\n");
        exit(1);
    }
    if(N < 1){
        printf("Error: N must be a positive integer...\n");
        exit(1);
    }
    if(L <= 0){
        printf("Error: L must be greater than 0...\n");
        exit(1);
    }

    printf("Allocating memory...\n");
    fgn = (double *) calloc(N*2, sizeof(double));
    fgn2 = (double *) calloc(N*2, sizeof(double));
    row_component = (double *) calloc(N, sizeof(double));
    reverse_component = (double *) calloc(N, sizeof(double));
    row = (double *) calloc(N*2, sizeof(double));
    outreal = (double *) calloc(N*2, sizeof(double));
    outimag = (double *) calloc(N*2, sizeof(double));
    w = (double *) calloc(N*2, sizeof(double));

    if(fgn == NULL || row_component == NULL || reverse_component == NULL || row == NULL){
        printf("Error: Allocating memory\n");
        exit(1);
    }

    srand(time(NULL));

    printf("Computing a normal distribution...\n");
    for(i=0; i<N; i++){
        rnd = (double)rand() / (double)RAND_MAX ;
        fgn[i] = 0.0 + 1.0*randomNormal();
    }

    if(H == 0.5){
        printf("Normal distribution has a H of 0.5, returning it...\n");
        printf("Writing to file...\n");
        for(i=0; i<N; i++){;
            fprintf(FP, "%lf\n", fgn[i]);
        }
    }
    else{
        j = N-2;
        for(i=0; i<N; i++){
            row_component[i] = autocovariance(H, i+1);
            reverse_component[j] = row_component[i];
            j--;
        }

        row[0] = 1.0;
        for(i=1; i<=N-1; i++){
            row[i] = row_component[i-1];
        }

        j = 0;
        row[N] = 0.0;
        for(i=N+1; i<N*2; i++){
            row[i] = reverse_component[j];
            j++;
        }

        dft(N*2, row, outreal, outimag, N*2);

        for(i=0; i<N*2; i++){
            rnd = (double)rand() / (double)RAND_MAX;
            fgn2[i] = 0.0 + 1.0 * randomNormal();
        }

        for(i=0; i<N*2; i++){
            if(outreal[i] < 0){
                printf("Error: Can't compute sequence for high H values and small length\n");
                flag = 0;
                break;
            }

            if(i == 0){
                w[i] = sqrt((outreal[i]) / (2.0 * N)) * fgn[i];
            }
            else if(i < N){
                w[i] = sqrt((outreal[i]) / (4.0 * N)) * (fgn[i] + 1 * fgn2[i]);
            }
            else if(i == N){
                w[i] = sqrt((outreal[i]) / (2.0 * N)) * fgn2[0];
            }
            else{
                w[i] = sqrt((outreal[i]) / (4.0 * N)) * (fgn[2 * N - i] - 1 * fgn2[2 * N - i]);
            }
        }

        if(flag == 1){
            dft(N*2, w, outreal, outimag, N*2);

            printf("Writing to file...\n", H, N);
            for(i=0; i<N; i++){
                fgn[i] = outreal[i];
                fgn[i] *= scale;
                fprintf(FP, "%f\n", fgn[i]);
            }
        }
    }

    printf("Freeing memory...\n");
    fclose(FP);
    free(fgn);
    free(fgn2);
    free(row_component);
    free(reverse_component);
    free(row);
    free(w);
    printf("Ending...\n");

    return 0;
}
