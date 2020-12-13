#include <utility>
#include <iostream>
#include <vector>
#include <assert.h>
#include <math.h>
#include <time.h>

using namespace std;

class LinearInterpolation {
public:
    LinearInterpolation(const vector<double> &t_via, const vector<double> &p_via) {
        if(p_via.size() != t_via.size()) throw runtime_error("The p_via and t_via must have a same length");
        p_via_ = p_via;
        t_via_ = t_via;
    }

    static pair<double, double> linear(double q0, double q1, double t0, double t1) {
        if(fabs(t0 - t1) < 1e-6) throw runtime_error("t0 and t1 must be different");
        double a0 = q0;
        double a1 = (q1 - q0)/(t1 - t0);
        return {a0, a1};
    }

    vector<double> getPosition(double t) {
        if((t < t_via_[0]) || (t > t_via_[t_via_.size() - 1])) throw runtime_error("The specific time error, time ranges error");
        // find the index of t1
        int i = 0, j = 0;
        for (; i < t_via_.size(); i++) {
            if (t_via_[i] >= t) {
                j = i;
                break;
            }
        }

        if(j == 0) {
            i = 0;
            j = 1;
        } else {
            i = j - 1;
        }
        
        vector<double> q(3,0);

        // position
        double q0 = p_via_[i];
        double t0 = t_via_[i];
        double q1 = p_via_[j];
        double t1 = t_via_[j];

        pair<double, double> ret = linear(q0, q1, t0, t1);
        double a0 = ret.first, a1 = ret.second;
        q[0] = a0 + a1*(t - t0);

        // velocity
        q[1] = a1;

        // acceleration
        q[2] = 0; // for linear model, the acceleration is infinite, here we set to zero
        return q;
    }

private:
    vector<double> p_via_;
    vector<double> t_via_;
};

class ParabolicInterpolation {
public:
    ParabolicInterpolation(const vector<double> &t_via, const vector<double> &p_via, const vector<double> &v_via) {
        if(p_via.size() != t_via.size()) throw runtime_error("The p_via and t_via must have a same length");
        if(v_via.size() != t_via.size()) throw runtime_error("The v_via and t_via must have a same length");
        p_via_ = p_via;
        v_via_ = v_via;
        t_via_ = t_via;
    }

    static vector<double> parabolic(double q0, double q1, double v0, double v1, double t0, double t1, double tf, double qf) {
        if(fabs(t0 - t1) < 1e-6) throw runtime_error("t0 and t1 must be different");
        if((tf <= t0) or (tf >= t1)) throw runtime_error("tf must satisfy t0 < tf < t1");
        if((qf <= min(q0, q1)) or (qf >= max(q0, q1))) throw runtime_error("qf must satisfy min(q0, q1) < qf < max(q0, q1)");

        double T = t1 - t0;
        double h = q1 - q0;
        double Ta = tf - t0;
        double Td = t1 - tf;

        double a0 = q0;
        double a1 = v0;
        double a2 = (2*h - v0*(T + Ta) - v1*Td)/(2*T*Ta);
        double a3 = (2*q1*Ta + Td*(2*q0 + Ta*(v0 - v1)))/(2*T);
        double a4 = (2*h - v0*Ta - v1*Td)/T;
        double a5 = -(2*h - v0*Ta - v1*(T+Td))/(2*T*Td);

        return {a0, a1, a2, a3, a4, a5};
    }

    vector<double> getPosition(double t) {
        if((t < t_via_[0]) || (t > t_via_[t_via_.size() - 1])) throw runtime_error("The specific time error, time ranges error");
        // find the index of t1
        int i = 0, j = 0;
        for (; i < t_via_.size(); i++) {
            if (t_via_[i] >= t) {
                j = i;
                break;
            }
        }

        if(j == 0) {
            i = 0;
            j = 1;
        } else {
            i = j - 1;
        }

        vector<double> q(3,0);

        // get given position
        double q0 = p_via_[i];
        double v0 = v_via_[i];
        double t0 = t_via_[i];

        double q1 = p_via_[j];
        double v1 = v_via_[j];
        double t1 = t_via_[j];

        // symmetric acceleration
        double tf = (t0 + t1)/2;
        double qf = (q0 + q1)/2;

        // asymmetric acceleration, specify tf and qf by users
        // tf = ?
        // qf = ?

        vector<double> ret = parabolic(q0, q1, v0, v1, t0, t1, tf, qf);
        double a0 = ret[0], a1 = ret[1], a2 = ret[2], a3 = ret[3], a4 = ret[4], a5 = ret[5];

        if(t <= tf) {
            q[0] = a0 + a1 * (t - t0) + a2 * pow((t - t0), 2);
            q[1] = a1 + 2 * a2 * (t - t0);
            q[2] = 2 * a2;
        } else {
            q[0] = a3 + a4 * (t - tf) + a5 * pow((t - tf), 2);
            q[1] = a4 + 2 * a5 * (t - tf);
            q[2] = 2 * a5;
        }

        return q;
    }

private:
    vector<double> p_via_;
    vector<double> v_via_;
    vector<double> t_via_;
};

class CubicInterpolation {
public:
    CubicInterpolation(const vector<double> &t_via, const vector<double> &p_via, const vector<double> &v_via) {
        if(p_via.size() != t_via.size()) throw runtime_error("The p_via and t_via must have a same length");
        if(v_via.size() != t_via.size()) throw runtime_error("The v_via and t_via must have a same length");
        p_via_ = p_via;
        v_via_ = v_via;
        t_via_ = t_via;
    }

    static vector<double> cubic(double q0, double q1, double v0, double v1, double t0, double t1) {
        if(fabs(t0 - t1) < 1e-6) throw runtime_error("t0 and t1 must be different");

        double T = t1 - t0;
        double h = q1 - q0;

        double a0 = q0;
        double a1 = v0;
        double a2 = (3*h - (2*v0 + v1)*T) / (pow(T, 2));
        double a3 = (-2*h + (v0 + v1)*T) / (pow(T, 3));

        return {a0, a1, a2, a3};
    }

    vector<double> getPosition(double t) {
        if((t < t_via_[0]) || (t > t_via_[t_via_.size() - 1])) throw runtime_error("The specific time error, time ranges error");
        // find the index of t1
        int i = 0, j = 0;
        for (; i < t_via_.size(); i++) {
            if (t_via_[i] >= t) {
                j = i;
                break;
            }
        }

        if(j == 0) {
            i = 0;
            j = 1;
        } else {
            i = j - 1;
        }

        vector<double> q(3,0);

        // get given position
        double q0 = p_via_[i];
        double v0 = v_via_[i];
        double t0 = t_via_[i];

        double q1 = p_via_[j];
        double v1 = v_via_[j];
        double t1 = t_via_[j];

        vector<double> ret = cubic(q0, q1, v0, v1, t0, t1);
        double a0 = ret[0], a1 = ret[1], a2 = ret[2], a3 = ret[3];

        q[0] = a0 + a1*(t - t0) + a2*pow((t - t0), 2) + a3*pow((t - t0), 3); // position
        q[1] = a1 + 2*a2*(t - t0) + 3*a3*pow((t - t0), 2); // velocity
        q[2] = 2*a2 + 6*a3*(t - t0); // acceleration

        return q;
    }

private:
    vector<double> p_via_;
    vector<double> v_via_;
    vector<double> t_via_;
};

class Polynomial5Interpolation {
public:
    Polynomial5Interpolation(const vector<double> &t_via, const vector<double> &p_via, const vector<double> &v_via, const vector<double> &a_via) {
        if(p_via.size() != t_via.size()) throw runtime_error("The p_via and t_via must have a same length");
        if(v_via.size() != t_via.size()) throw runtime_error("The v_via and t_via must have a same length");
        if(a_via.size() != t_via.size()) throw runtime_error("The a_via and t_via must have a same length");
        p_via_ = p_via;
        v_via_ = v_via;
        a_via_ = a_via;
        t_via_ = t_via;
    }

    static vector<double> polynomial(double q0, double q1, double v0, double v1, double acc0, double acc1, double t0, double t1) {
        if(fabs(t0 - t1) < 1e-6) throw runtime_error("t0 and t1 must be different");

        double T = t1 - t0;
        double h = q1 - q0;

        double a0 = q0;
        double a1 = v0;
        double a2 = acc0/2;
        double a3 = (20*h - (8*v1 + 12*v0)*T - (3*acc0 - acc1)*pow(T, 2)) / (2*pow(T, 3));
        double a4 = (-30*h + (14*v1 + 16*v0)*T + (3*acc0 - 2*acc1)*pow(T, 2)) / (2*pow(T, 4));
        double a5 = (12*h - 6*(v1 + v0)*T + (acc1 - acc0)*pow(T, 2)) / (2*pow(T, 5));

        return {a0, a1, a2, a3, a4, a5};
    }

    vector<double> getPosition(double t) {
        if((t < t_via_[0]) || (t > t_via_[t_via_.size() - 1])) throw runtime_error("The specific time error, time ranges error");
        // find the index of t1
        int i = 0, j = 0;
        for (; i < t_via_.size(); i++) {
            if (t_via_[i] >= t) {
                j = i;
                break;
            }
        }

        if(j == 0) {
            i = 0;
            j = 1;
        } else {
            i = j - 1;
        }

        vector<double> q(3,0);

        // get given position
        double q0 = p_via_[i];
        double v0 = v_via_[i];
        double acc0 = a_via_[i];
        double t0 = t_via_[i];

        double q1 = p_via_[j];
        double v1 = v_via_[j];
        double acc1 = a_via_[j];
        double t1 = t_via_[j];

        vector<double> ret = polynomial(q0, q1, v0, v1, acc0, acc1, t0, t1);
        double a0 = ret[0], a1 = ret[1], a2 = ret[2], a3 = ret[3], a4 = ret[4], a5 = ret[5];

        q[0] = a0 + a1*(t - t0) + a2*pow((t - t0), 2) + a3*pow((t - t0), 3) + a4*pow((t - t0), 4) + a5*pow((t - t0), 5); // position
        q[1] = a1 + 2*a2*(t - t0) + 3*a3*pow((t - t0), 2) + 4*a4*pow((t - t0), 3) + 5*a5*pow((t - t0), 4); // velocity
        q[2] = 2*a2 + 6*a3*(t - t0) + 12*a4*pow((t - t0), 2) + 20*a5*pow((t - t0), 3); // acceleration

        return q;
    }

private:
    vector<double> p_via_;
    vector<double> v_via_;
    vector<double> a_via_;
    vector<double> t_via_;
};



int main() {
    vector<double> t_given = {0, 1, 3, 4.5, 6, 8, 10};
    vector<double> p_given = {0, 1.6, 3.2, 2, 4, 0.2, 1.2};
    vector<double> v_given = {0, 1.0, 2.0, -2.0, -1.0, 0, 0};
    vector<double> a_given = {0, 1, 2, 3, 2, 1, 0};

    LinearInterpolation linear_interpolation(t_given, p_given);
    ParabolicInterpolation parabolic_interpolation(t_given, p_given, v_given);
    CubicInterpolation cubic_interpolation(t_given, p_given, v_given);
    Polynomial5Interpolation polynomial5_interpolation(t_given, p_given, v_given, a_given);

    const int cnt = 1000;
    vector<double> t(cnt);
    double step = (t_given[t_given.size() - 1] - t_given[0]) / cnt;
    for(int i = 0; i < cnt; i++)
        t[i] = t_given[0] + step*i;

    clock_t start = clock();
    for(auto i : t) {
//        printf("%f\n", linear_interpolation.getPosition(ti)[0]);
//        printf("%f\n", parabolic_interpolation.getPosition(ti)[0]);
//        printf("%f\n", cubic_interpolation.getPosition(ti)[2]);
//        printf("%f\n", polynomial5_interpolation.getPosition(i)[2]);


        linear_interpolation.getPosition(i); // 0.179ms
//        parabolic_interpolation.getPosition(i); // 0.329ms
//        cubic_interpolation.getPosition(i); // 0.447ms
//        polynomial5_interpolation.getPosition(i); // 0.839ms
    }

    cout << "time[1000 points] = "<< double(clock()-start) / CLOCKS_PER_SEC * 1000 << "ms" << endl;

    return 0;
}