#include <iostream>
#include <vector>
#include <bits/stdc++.h>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#include "stdio.h"
#include "x86intrin.h"
#include <limits>
#include 	"x86intrin.h"
#include 	"stdio.h"

#define		VECTOR_SIZE		  	4194304

using namespace std;

float *generateArr() {
    float *arr;
	arr = new float [VECTOR_SIZE];

	for (long i = 0; i < VECTOR_SIZE; i++)
		arr[i] = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/100.0));
    return arr;
}

float serialQ1(float *arr){
    struct timeval start, end;
	gettimeofday(&start, NULL);
    
    float min = arr[0];
    int indx = 0;

    for (long i=1; i<VECTOR_SIZE;i++){
        if (arr[i] < min){
            indx = i;
            min = arr[i];
        }
    }

    gettimeofday(&end, NULL);
    long second = end.tv_sec - start.tv_sec;
    long micros = ((second * 1000000) + end.tv_usec) - start.tv_usec;

    cout << "Minimum of array in serial execution is " << min << " and the index is " << indx << endl;
    cout << "Serial execution time is " << second << " second and " << micros << " micro second" << endl << endl;
    return (second * 1000000) + micros;

}

union float_vec {
    __m128 v;
    float vec[4];
};

union int_vec {
    __m128i v;
    int vec[4];
};

float parallelQ1(float *arr){
    struct timeval start, end;
	gettimeofday(&start, NULL);

    __m128 vec, smaller;
    __m128i indexes = _mm_setr_epi32(0, 1, 2, 3);
    __m128i Incr = _mm_setr_epi32(4, 4, 4, 4);
    __m128i min_indexes = indexes;
    __m128 minimumValues = _mm_load_ps(arr);
    __m128i iSmaller;

    for (int i = 0; i < VECTOR_SIZE; i+=4) {
        vec = _mm_loadu_ps(&arr[i]);
        smaller = _mm_cmplt_ps(vec, minimumValues);
        iSmaller = __m128i(smaller);
        min_indexes = _mm_blendv_epi8(min_indexes, indexes, iSmaller);
        minimumValues = _mm_min_ps(vec, minimumValues);
        indexes = _mm_add_epi32(indexes, Incr);
    }

    float_vec f_vec = {minimumValues};

    int minIndx = 0;
    for (int i = 1; i < 4; i++)
        if (f_vec.vec[i] < f_vec.vec[minIndx])
            minIndx = i;

    int_vec idx_vec = {min_indexes};
    int finalMinIndx = int(idx_vec.vec[minIndx]);

    gettimeofday(&end, NULL);
    long second = end.tv_sec - start.tv_sec;
    long micros = ((second * 1000000) + end.tv_usec) - start.tv_usec;

    cout << "Minimum of array in parallel execution is " << f_vec.vec[minIndx] << " and index is " << finalMinIndx << endl;
    cout << "Parallel execution time is " << second << " second and " << micros << " micro second" << endl << endl;

    return (second * 1000000) + micros;
}

float serialQ2(float *arr){
    struct timeval start, end;
	gettimeofday(&start, NULL);
    
    float sum, avg, std;

    float tempStd[4];
    float tempResult[4];

    tempResult[0] = tempResult[1] = tempResult[2] = tempResult[3] = 0.0;
    tempStd[0] = tempStd[1] = tempStd[2] = tempStd[3] = 0.0;

	for (long i = 0; i < VECTOR_SIZE; i+=4)
		tempResult[0] += arr[i];
	for (long i = 0; i < VECTOR_SIZE; i+=4)
		tempResult[1] += arr[i + 1];
	for (long i = 0; i < VECTOR_SIZE; i+=4)
		tempResult[2] += arr[i + 2];
	for (long i = 0; i < VECTOR_SIZE; i+=4)
		tempResult[3] += arr[i + 3];

	sum = tempResult[0] + tempResult[1] + tempResult[2] + tempResult[3];

    avg = sum / VECTOR_SIZE;

    for (long i = 0; i < VECTOR_SIZE; i+=4)
		tempStd[0] += (arr[i] - avg) * (arr[i] - avg);
	for (long i = 0; i < VECTOR_SIZE; i+=4)
		tempStd[1] += (arr[i + 1] - avg) * (arr[i + 1] - avg);
	for (long i = 0; i < VECTOR_SIZE; i+=4)
		tempStd[2] += (arr[i + 2] - avg) * (arr[i + 2] - avg);
	for (long i = 0; i < VECTOR_SIZE; i+=4)
		tempStd[3] += (arr[i + 3] - avg) * (arr[i + 3] - avg);

    std = tempStd[0] + tempStd[1] + tempStd[2] + tempStd[3];

    std = sqrt(std/VECTOR_SIZE);

    gettimeofday(&end, NULL);
    long second = end.tv_sec - start.tv_sec;
    long micros = ((second * 1000000) + end.tv_usec) - start.tv_usec;

    cout << "Average of array in serial execution is " << avg << " and std is " << std << endl;
    cout << "Serial execution time is " << second << " second and " << micros << " micro second" << endl << endl;
    return (second * 1000000) + micros;
}

__m128 averageParallel(float *arr){
	__m128 vec;
	__m128 sum = _mm_set1_ps(0.0f);
    __m128 avg = _mm_set1_ps(VECTOR_SIZE);
	for (long i = 0; i < VECTOR_SIZE; i+=4) {
		vec = _mm_loadu_ps(&arr[i]);
		sum = _mm_add_ps(sum, vec);
	}
	sum = _mm_hadd_ps(sum, sum);
	sum = _mm_hadd_ps (sum, sum);

    avg = _mm_div_ps(sum, avg);

    return avg;
}


__m128 stdParallel(float *arr, __m128 avg){
	__m128 vec, sub, pow, std, sum;

    std = _mm_set1_ps(VECTOR_SIZE);
    sum = _mm_set1_ps(0.0f);
	for (long i = 0; i < VECTOR_SIZE; i+=4) {
		vec = _mm_loadu_ps(&arr[i]);
		sub = _mm_sub_ps(vec, avg);
        pow = _mm_mul_ps(sub,sub);
		sum = _mm_add_ps(sum, pow);
	}
    
	sum = _mm_hadd_ps(sum, sum);
	sum = _mm_hadd_ps (sum, sum);

    std = _mm_div_ps(sum, std);
    std = _mm_sqrt_ps(std);
    return std;
}


float parallelQ2(float *arr){
    struct timeval start, end;
	gettimeofday(&start, NULL);

    __m128 avg = averageParallel(arr);
    float avgFloat = _mm_cvtss_f32 (avg);

    __m128 std = stdParallel(arr, avg);
    float stdFloat = _mm_cvtss_f32 (std);

    gettimeofday(&end, NULL);
    long second = end.tv_sec - start.tv_sec;
    long micros = ((second * 1000000) + end.tv_usec) - start.tv_usec;

    cout << "Average of array in parallel execution is " << avgFloat << " and std is " << stdFloat << endl;
    cout << "Parallel execution time is " << second << " second and " << micros << " micro second" << endl << endl;

    return (second * 1000000) + micros;
}



int main(){
    float *arr = generateArr();
    float serialTime, parallelTime;
    cout << endl << "Group Members" << endl;
    cout << "Kiavash Jamshidi: " << 810197486 << endl;
    cout << "Sina Salimian: " << 810197527 << endl << endl;

    serialTime = serialQ1(arr);
    parallelTime = parallelQ1(arr);
    cout << "Speedup for Question1 is " << serialTime / parallelTime << endl << endl;

    serialTime = serialQ2(arr);
    parallelTime = parallelQ2(arr);
    cout << "Speedup for Question2 is " << serialTime / parallelTime << endl << endl;
}