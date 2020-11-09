#ifndef _TRAJECTORY_GENERATOR_H_
#define _TRAJECTORY_GENERATOR_H_

#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>
// #include "mosek.h"
#include "osqp.h"
#include "bezier_base.h"
#include "data_type.h"

using namespace std;
using namespace Eigen;

class TrajectoryGenerator {
private:

public:
        TrajectoryGenerator(){}
        ~TrajectoryGenerator(){}

        /* Use Bezier curve for the trajectory */
       int BezierPloyCoeffGeneration(
            const vector<Cube> &corridor,
            const MatrixXd &MQM,
            const MatrixXd &pos,
            const MatrixXd &vel,
            const MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const int traj_order,
            const double minimize_order,
            const double margin,
            const bool & isLimitVel,
            const bool & isLimitAcc,
            double & obj,
            MatrixXd & PolyCoeff);
protected:

        virtual OSQPSettings* SolverDefaultSettings();
        
        void FreeData(OSQPData* data);

        template <typename T> T* CopyData(const std::vector<T>& vec) {
        T* data = new T[vec.size()];
        memcpy(data, vec.data(), sizeof(T) * vec.size());
        return data;
    }
};

#endif
