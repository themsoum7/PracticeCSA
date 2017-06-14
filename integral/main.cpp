#include<iostream>
#include <cmath>
#include <thread>
#include <mutex>
#include <fstream>
#include <chrono>
#include <atomic>
#include <math.h>

using namespace std;

double temp;
double res1;
double res2 = 0;
double a = 20;
double b = 0.2;
double c = 2 * M_PI;
mutex mtx;

float integral_func(float x1, float x2){
    double res;
    res = -a * exp(-b * sqrt(1/2 * ((pow(x1, 2) + pow(x2, 2))))) - exp(1/2*(cos(c*x1) + cos(c*x2))) + a + exp(1);
    return (float) res;
}

//Tried something another than the final version.

/*double sumIntegral(double lowbound, int n, double dx){
    double cumSum = 0;
    for(int i = 0; i < n; i++){
        double xi = lowbound + i * dx;
        double funValue = integral_func(xi, );
        double area = funValue * dx;
        cumSum += area;
    }
    return cumSum;
}

double sumDoubleIntegral(double lowbound1, double lowbound2, int n, double dy){
    double cumSum2 = 0;
    for(int i = 0; i < n; i++){
        double yi = lowbound1 + i * dy;
        double dx = (yi - lowbound2) / n;
        double smallSum = 0;

        for(int j = 0; j < n; j++){
            double xi = lowbound2 + dx * j;
            double funValue = ;
            double area = funValue * dx;
            smallSum += area;
        }
        double secondArea = smallSum * dy;
        cumSum2 += secondArea;
    }
    return cumSum2;
}*/

double calc_integral(double lowbound1, double lowbound2, int m, int l, double dx, double dy){
    double Sum = 0;
    double xi = lowbound1;
    for(int i = 0; i < m; i++){
        double yi = lowbound2;
        xi += dx;
        for(int j = 0; j < l; j++){
            double xv = (2 * xi + dx) / 2;
            double yv = (2 * yi + dy) / 2;
            double func_res = integral_func((float)xv, (float)yv);
            func_res *= dx * dy;
            Sum += func_res;
            yi += dy;
        }
    }

    lock_guard<mutex> lguard(mtx);
    res1 += Sum;

    if (res2 == 0){
        temp = Sum;
        res2 = temp;
        res1 = 0;
    }
}

inline chrono::high_resolution_clock::time_point get_current_time_fenced() {
    atomic_thread_fence(memory_order_seq_cst);
    auto res_time = chrono::high_resolution_clock::now();
    atomic_thread_fence(memory_order_seq_cst);
    return res_time;
}

template<class D>
inline long long to_us(const D& d) {
    return chrono::duration_cast<chrono::microseconds>(d).count();
}

int main(int a0rgc, char *argv[]) {
    double lowerbound_exp;
    double upperbound_exp;
    double lowerbound_exp2;
    double upperbound_exp2;
    cout<< "x-lower bound: ";
    cin >> lowerbound_exp;
    cout<< "x-upper bound: ";
    cin >> upperbound_exp;
    cout<< "y-lower bound: ";
    cin >> lowerbound_exp2;
    cout<< "y-upper bound: ";
    cin>>upperbound_exp2;

    cout << "intervals: ";
    int m;
    cin >> m;
    int l = m;
    double dx1 = (upperbound_exp - lowerbound_exp) / m;
    double dy1 = (upperbound_exp2 - lowerbound_exp2) / l;
    calc_integral(lowerbound_exp, lowerbound_exp2, m, l, dx1, dy1);

    cout<<"number of threads to use: ";
    int threads;
    cin >> threads;
    thread th[threads];

    double upperbound_a = upperbound_exp/threads;
    double lowerbound_temp = lowerbound_exp;
    double lowerbound_temp2 = lowerbound_exp2;
    double upperbound_temp = upperbound_a;
    double dx = (upperbound_a - lowerbound_exp) / m;
    double dy = (upperbound_exp2 - lowerbound_exp2) / l;


    auto begin = get_current_time_fenced();
    for (int i = 0; i < threads; ++i) {
        th[i] = thread(calc_integral, lowerbound_temp, lowerbound_temp2, m, l, dx, dy);
        lowerbound_temp += upperbound_a;
        upperbound_temp += upperbound_a;
    }
    for (int i = 0; i < threads; ++i) {
        th[i].join();
    }

    auto finish = get_current_time_fenced();
    auto time_used = finish - begin;

    ofstream f("file.txt");
    f << "Result of calculation: " << res1<< endl;
    f << "Expected error: " << res1 - res2 << endl;
    f << "Total time of calculation: " << to_us(time_used) << endl;
    return 0;
}