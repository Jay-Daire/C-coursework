//Including all the necessary extra libraries for this program, as well as defining PI. The first few times I ran this, the values returned were off,
//so I defined it to 50 figures.
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.14159265358979323846264338327950288419716939937510

//First I define complex numbers as a structure with two components, real and imaginary. I chose to use doubles as they have more space allocated to them, allowing
//for more accuracy.
typedef struct Complex {
    double re;
    double im;
} Complex;

//I then define h1 and h2 as functions that can be generated from values of t.
Complex h1(double t) {
    Complex h1out;
    h1out.re= cos(t) + cos(5.0 * t);
    h1out.im= sin(t) + sin(5.0 * t);
    return h1out;
}
Complex h2(double t) {
    Complex h2out;
    h2out.re= exp(pow((t - PI), 2.0) / 2.0);
    h2out.im= 0;
    return h2out;
}

//I designed a function similar to linspace in python, to create an array of evenly spaced values between two input values. Again, I use doubles for accuracy.
double* linspace(const double init, const double fin, const int N) {
    //Initialising all used variables.
    double spacing = (double) (fin - init) / N;
    double *arr = (double *) malloc(N * sizeof(double));
    double val_hold = init;
    //Looping from the initial value to final value N times, adding spacing value each time to produce the spaced values.
    for(int i = 0; i < N; ++i) {
        arr[i] = val_hold;
        val_hold += spacing;
    }
    return arr;
}

//This is my Discrete Fourier Transform function. it takes an array of complex values and runs it through the fourier transform,
//then outputs the transformed complex number array
Complex* dft(const Complex* h_in, const int N) {
    Complex* H_out = (Complex*)malloc(N * sizeof(Complex));
    Complex val_hold;
    //Here I use two for loops. The inner loop calculates the overall value for each space in the array I return, whereas the outer loop populates the array and resets
    //the placeholder value.
    for(int n = 0; n < N; ++n) {
        val_hold.re = 0;
        val_hold.im = 0;
        for(int k = 0; k < N; ++k) {
            double omega = -(2 * PI * n * k) / N;
            val_hold.re += (h_in[k].re * cos(omega)) - (h_in[k].im * sin(omega));
            val_hold.im += (h_in[k].im * cos(omega)) + (h_in[k].re * sin(omega));
        }
        H_out[n] = val_hold;
    }
    return H_out;
}

//This is the Inverse Fourier Transform function - with a built in skip function. It takes in an array of complex values and runs them through the inverse transform
//skipping the number input for skip.
Complex* ift(const Complex* H_in, const int N, const int skip) {
    Complex* h_inv_out = (Complex*)malloc(N * sizeof(Complex));
    Complex val_hold;
    //As with the DFT, the inner loop calculates the value to be put into the array, and the outer loop places the value into the array to be returned.
    //I also added an if condition so I could skip certain values as necessary.
    for(int k = 0; k < N; ++k) {
        val_hold.re = 0;
        val_hold.im = 0;
        for(int n = 0; n < N; ++n) {
            if(n == skip) {
                val_hold.re += 0;
                val_hold.im += 0;
            }
            else {
                double omega = (2 * PI * n * k) / N;
                val_hold.re += (H_in[n].re * cos(omega)) - (H_in[n].im * sin(omega));
                val_hold.im += (H_in[n].im * cos(omega)) + (H_in[n].re * sin(omega));
            }
        }
        h_inv_out[k].re = val_hold.re / N;
        h_inv_out[k].im = val_hold.im / N;
    }
    return h_inv_out;
}

//A second IFT function built specifically for H3. I had trouble figuring out how to generalise the original IFT function enough to use H3 in the way I needed to.
Complex* iftH_3(const Complex* H_in, const int N, const int inner_N, const int* index_in) {
    Complex* h_inv_out = (Complex*)malloc(N * sizeof(Complex));
    Complex val_hold;
    //The main difference from the original IFT is I removed the if condition, and pulled the index number of each H value from a sorted array I created.
    for(int k = 0; k < N; ++k) {
        val_hold.re = 0;
        val_hold.im = 0;
        for(int n = 0; n < inner_N; ++n) {
            double omega = (2 * PI * index_in[n] * k) / N;
            val_hold.re += (H_in[n].re * cos(omega)) - (H_in[n].im * sin(omega));
            val_hold.im += (H_in[n].im * cos(omega)) + (H_in[n].re * sin(omega));
        }
        h_inv_out[k].re = val_hold.re / N;
        h_inv_out[k].im = val_hold.im / N;
    }
    return h_inv_out;
}

//Created a read_function to tidy up the main program - this simply reads in a text file for as many rows as you ask it to and outputs a complex number array.
Complex* read_func(const char* filename, const int rows) {
    FILE* fp = fopen(filename, "r");
    double re, im;
    int i = 0;
    Complex* data = (Complex*)malloc(rows * sizeof(Complex));
    //In my while loop, I ignore the index and time columns, and have two conditions. i must be smaller than the amount of rows input, and I must not have reached the end
    //of the file.
    while (fscanf(fp, "%*d, %*lf, %lf, %lf\n", &re, &im) != EOF && i<rows) {
        data[i].re = re;
        data[i].im = im;
        ++i;
    }
    fclose(fp);
    return data;
}

//This is a quick function I built to find the modulus of a complex number. It takes in the complex number, squares the real and imaginary parts and sums them then
//find the square root of the sum, which is the modulus.
double modulus(const Complex input) {
    double mod_value;
    mod_value = sqrt(pow(input.re, 2.0) + pow(input.im, 2.0));
    return mod_value;
}

//This performs quicksort on an array of complex numbers, and records their initial positions
//for the complex quicksort I copied the built-in quicksort function for C, and added an index to be edited alongside the sorted H3 numbers.
int complex_quicksort(Complex* H_in, int* index, const int init, const int fin) {
    int i, j, pivot, index_holder;
    Complex H_holder;

    if (init < fin) {
        pivot = init;
        i = init;
        j = fin;

        while (i < j) {
            while (modulus(H_in[i]) <= modulus(H_in[pivot]) && i < fin)
                ++i;
            while (modulus(H_in[j]) > modulus(H_in[pivot]))
                --j;
            if (i < j) {
                H_holder = H_in[i];
                index_holder = index[i];
                H_in[i] = H_in[j];
                index[i] = index[j];
                H_in[j] = H_holder;
                index[j] = index_holder;
            }
        }
        H_holder = H_in[pivot];
        index_holder = index[pivot];
        H_in[pivot] = H_in[j];
        index[pivot] = index[j];
        H_in[j] = H_holder;
        index[j] = index_holder;
        complex_quicksort(H_in, index, init, j -1);
        complex_quicksort(H_in, index, j + 1, fin);
    }
    return 0;
}

int main() {
    //Friendly message to start off.
    printf("Hello, Dr Hendrik!\n");
    // Generating the times used for h1 and h2, then initialising an index for the complex quicksort later. I also generate the h3 time array.
    //As well as all that, I produce an H3 top 4 index values for my special IFT function.
    double* timespace = linspace(0, 2 * PI, 100);
    double* h3time = linspace(0, 2 * PI, 200);
    int* index_values = (int*)malloc(200);
    int* H_3_index4 = (int*)malloc(4);

    //Initialising an array for h1 and h2, then reading h3 from the file.
    Complex* h_1 = (Complex*)malloc(sizeof(Complex) * 100);
    Complex* h_2 = (Complex*)malloc(sizeof(Complex) * 100);
    Complex* h_3 = read_func("h3data.txt", 200);
    Complex* H_3_top4 = (Complex*)malloc(sizeof(Complex) * 4);

    //All of the file opening and writing to is done here - I tried to create a write function but it didn't work.
    FILE* fp1 = fopen("h1data.txt", "w");
    FILE* fp2 = fopen("h2data.txt", "w");
    FILE* fp3 = fopen("h'1data.txt", "w");
    FILE* fp4 = fopen("h'2data.txt", "w");
    FILE* fp5 = fopen("h'3data.txt", "w");
    fprintf(fp1,"time,re,im\n");
    fprintf(fp2,"time,re,im\n");
    fprintf(fp3,"time,re,im\n");
    fprintf(fp4,"time,re,im\n");
    fprintf(fp5, "time,re,im\n");

    //Filling the index values array with.. indexing values.
    for(int i = 0; i < 200; ++i) {
        index_values[i] = i;
    }

    //Here I'm writing the time and complex number for h1 and h2 to file.
    for(int i = 0; i < 100; ++i) {
        h_1[i] = h1(timespace[i]);
        h_2[i] = h2(timespace[i]);
        fprintf(fp1, "%lf,%lf,%lf\n", timespace[i], h_1[i].re, h_1[i].im);
        fprintf(fp2, "%lf,%lf,%lf\n", timespace[i], h_2[i].re, h_2[i].im);
    }

    //Now I create H1, H2 and H3. I had to do this down here as h_1, h_2 and h_3 don't exist properly until now. I also sort H_3 from smallest to largest values
    //and sort its corresponding index array, so I can do the IFT transform correctly.
    Complex* H_1 = dft(h_1, 100);
    Complex* H_2 = dft(h_2, 100);
    Complex* H_3 = dft(h_3, 200);
    complex_quicksort(H_3, index_values, 0, 200);
    for(int i = 0; i < 4; ++i) {
        H_3_top4[i] = H_3[199-i];
        H_3_index4[i] = index_values[199-i];
    }

    //Printing H1 and H2 values to console.
    printf("H1\n");
    printf("Real, Imaginary\n");
    for(int i = 0; i < 100; ++i) {
        printf("%lf,%lf\n", H_1[i].re, H_1[i].im);
    }
    printf("H2\n");
    printf("Real, Imaginary\n");
    for(int i = 0; i < 100; ++i) {
        printf("%lf,%lf\n", H_2[i].re, H_2[i].im);
    }

    //Producing the h'1, h'2 and h'3 values using my IFT functions.
    Complex* h_inv1 = ift(H_1, 100, 1);
    Complex* h_inv2 = ift(H_2, 100, 0);
    Complex* h_inv3 = iftH_3(H_3_top4, 200, 4, H_3_index4);

    //Here I write all the inverse transform results to file so I can plot them.
    for(int i = 0; i < 100; ++i) {
        fprintf(fp3, "%lf,%lf,%lf\n", timespace[i], h_inv1[i].re, h_inv1[i].im);
        fprintf(fp4, "%lf,%lf,%lf\n", timespace[i], h_inv2[i].re, h_inv2[i].im);
    }
    for(int i = 0; i < 200; ++i) {
        fprintf(fp5, "%lf,%lf,%lf\n", h3time[i], h_inv3[i].re, h_inv3[i].im);
    }

    //Closing all the files for good practice
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fp4);
    fclose(fp5);
    return 0;
}
