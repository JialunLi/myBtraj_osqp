#include "trajectory_generator.h"
using namespace std;
using namespace Eigen;

int TrajectoryGenerator::BezierPloyCoeffGeneration(
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
            MatrixXd & PolyCoeff)  // define the order to which we minimize.   1 -- velocity, 2 -- acceleration, 3 -- jerk, 4 -- snap  
{   
#define ENFORCE_VEL  isLimitVel // whether or not adding extra constraints for ensuring the velocity feasibility
#define ENFORCE_ACC  isLimitAcc // whether or not adding extra constraints for ensuring the acceleration feasibility

    double initScale = corridor.front().t;
    double lstScale  = corridor.back().t;
    int segment_num  = corridor.size();

    int n_poly = traj_order + 1;
    int s1d1CtrlP_num = n_poly;
    int s1CtrlP_num   = 3 * s1d1CtrlP_num; //number of coef each segment

    int equ_con_s_num = 3 * 3; // p, v, a in x, y, z axis at the start point
    int equ_con_e_num = 3 * 3; // p, v, a in x, y, z axis at the end point
    int equ_con_continuity_num = 3 * 3 * (segment_num - 1);
    int equ_con_num   = equ_con_s_num + equ_con_e_num + equ_con_continuity_num; // p, v, a in x, y, z axis in each segment's joint position
    
    int vel_con_num = 3 *  traj_order * segment_num;
    int acc_con_num = 3 * (traj_order - 1) * segment_num;

    if( !ENFORCE_VEL )
        vel_con_num = 0;

    if( !ENFORCE_ACC )
        acc_con_num = 0;

    int high_order_con_num = vel_con_num + acc_con_num; 
   
    int ctrlP_num = segment_num * s1CtrlP_num;

    vector<c_float> A_data;
    vector<c_int> A_indices;
    vector<c_int> A_indptr;
    vector<c_float> lower_bounds;
    vector<c_float> upper_bounds;
    
    /*** constraints bound ***/
    if(ENFORCE_VEL)
    {
        /***  Stack the bounding value for the linear inequality for the velocity constraints  ***/
        for(int i = 0; i < vel_con_num; i++)
        {
            lower_bounds.push_back(- maxVel);
            upper_bounds.push_back(+ maxVel);
        }
    }

    if(ENFORCE_ACC)
    {
        /***  Stack the bounding value for the linear inequality for the acceleration constraints  ***/
        for(int i = 0; i < acc_con_num; i++)
        {
            lower_bounds.push_back(- maxAcc);
            upper_bounds.push_back(+ maxAcc);
        }
    }

    //ROS_WARN("[Bezier Trajectory] equality bound %d", equ_con_num);
    for(int i = 0; i < equ_con_num; i ++ ){ 
        double beq_i;
        if(i < 3)                    beq_i = pos(0, i);
        else if (i >= 3  && i < 6  ) beq_i = vel(0, i - 3); 
        else if (i >= 6  && i < 9  ) beq_i = acc(0, i - 6);
        else if (i >= 9  && i < 12 ) beq_i = pos(1, i - 9 );
        else if (i >= 12 && i < 15 ) beq_i = vel(1, i - 12);
        else if (i >= 15 && i < 18 ) beq_i = acc(1, i - 15);
        else beq_i = 0.0;

        lower_bounds.push_back(beq_i);
        upper_bounds.push_back(beq_i);
    }

    
    //ROS_WARN("[Bezier Trajectory] Start stacking the Linear Matrix A, inequality part");
    vector< vector < pair <c_int,c_float> >> variables(ctrlP_num);
    int constraint_index = 0;
    // The velocity constraints
    if(ENFORCE_VEL)
    {   
        for(int k = 0; k < segment_num ; k ++ )
        {   
            for(int i = 0; i < 3; i++)
            {  // for x, y, z loop
                for(int p = 0; p < traj_order; p++)
                {
                    int nzi = 2;
                    int asub[nzi];
                    double aval[nzi];

                    aval[0] = -1.0 * traj_order;
                    aval[1] =  1.0 * traj_order;

                    asub[0] = k * s1CtrlP_num + i * s1d1CtrlP_num + p;    
                    asub[1] = k * s1CtrlP_num + i * s1d1CtrlP_num + p + 1;

                    variables[asub[0]].emplace_back(constraint_index, aval[0]);
                    variables[asub[1]].emplace_back(constraint_index, aval[1]);
                    constraint_index ++;
                }
            }
        }
    }

    // The acceleration constraints
    if(ENFORCE_ACC)
    {
        for(int k = 0; k < segment_num ; k ++ )
        {
            for(int i = 0; i < 3; i++)
            { 
                for(int p = 0; p < traj_order - 1; p++)
                {    
                    int nzi = 3;
                    int asub[nzi];
                    double aval[nzi];

                    aval[0] =  1.0 * traj_order * (traj_order - 1) / corridor[k].t;
                    aval[1] = -2.0 * traj_order * (traj_order - 1) / corridor[k].t;
                    aval[2] =  1.0 * traj_order * (traj_order - 1) / corridor[k].t;
                    asub[0] = k * s1CtrlP_num + i * s1d1CtrlP_num + p;    
                    asub[1] = k * s1CtrlP_num + i * s1d1CtrlP_num + p + 1;    
                    asub[2] = k * s1CtrlP_num + i * s1d1CtrlP_num + p + 2;    
                    
                    variables[asub[0]].emplace_back(constraint_index, aval[0]);
                    variables[asub[1]].emplace_back(constraint_index, aval[1]);
                    variables[asub[2]].emplace_back(constraint_index, aval[2]);
                    constraint_index ++;
                }
            }
        }
    }
    /*   Start position  */
    {
        // position :
        for(int i = 0; i < 3; i++)
        {  // loop for x, y, z       
            int nzi = 1;
            int asub[nzi];
            double aval[nzi];
            aval[0] = 1.0 * initScale;
            asub[0] = i * s1d1CtrlP_num;

            variables[asub[0]].emplace_back(constraint_index, aval[0]);
            
            constraint_index ++;
        }
        // velocity :
        for(int i = 0; i < 3; i++)
        {  // loop for x, y, z       
            int nzi = 2;
            int asub[nzi];
            double aval[nzi];
            aval[0] = - 1.0 * traj_order;
            aval[1] =   1.0 * traj_order;
            asub[0] = i * s1d1CtrlP_num;
            asub[1] = i * s1d1CtrlP_num + 1;

            variables[asub[0]].emplace_back(constraint_index, aval[0]);
            variables[asub[1]].emplace_back(constraint_index, aval[1]);
            
            constraint_index ++;
        }
        // acceleration : 
        for(int i = 0; i < 3; i++)
        {  // loop for x, y, z       
            int nzi = 3;
            int asub[nzi];
            double aval[nzi];
            aval[0] =   1.0 * traj_order * (traj_order - 1) / initScale;
            aval[1] = - 2.0 * traj_order * (traj_order - 1) / initScale;
            aval[2] =   1.0 * traj_order * (traj_order - 1) / initScale;
            asub[0] = i * s1d1CtrlP_num;
            asub[1] = i * s1d1CtrlP_num + 1;
            asub[2] = i * s1d1CtrlP_num + 2;

            variables[asub[0]].emplace_back(constraint_index, aval[0]);
            variables[asub[1]].emplace_back(constraint_index, aval[1]);
            variables[asub[2]].emplace_back(constraint_index, aval[2]);
            //r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            //row_idx ++;
            constraint_index ++;
        }
    }      

    /*   End position  */
    //ROS_WARN(" end position");
    {   
        // position :
        for(int i = 0; i < 3; i++)
        {  // loop for x, y, z       
            int nzi = 1;
            int asub[nzi];
            double aval[nzi];
            asub[0] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num;
            aval[0] = 1.0 * lstScale;
            // r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            // row_idx ++;
            variables[asub[0]].emplace_back(constraint_index, aval[0]);
            constraint_index ++;
        }
        // velocity :
        for(int i = 0; i < 3; i++)
        { 
            int nzi = 2;
            int asub[nzi];
            double aval[nzi];
            asub[0] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num - 1;
            asub[1] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num;
            aval[0] = - 1.0;
            aval[1] =   1.0;
            // r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            // row_idx ++;
            variables[asub[0]].emplace_back(constraint_index, aval[0]);
            variables[asub[1]].emplace_back(constraint_index, aval[1]);
            constraint_index ++;
        }
        // acceleration : 
        for(int i = 0; i < 3; i++)
        { 
            int nzi = 3;
            int asub[nzi];
            double aval[nzi];
            asub[0] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num - 2;
            asub[1] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num - 1;
            asub[2] = ctrlP_num - 1 - (2 - i) * s1d1CtrlP_num;
            aval[0] =   1.0 / lstScale;
            aval[1] = - 2.0 / lstScale;
            aval[2] =   1.0 / lstScale;
            // r = MSK_putarow(task, row_idx, nzi, asub, aval);    
            // row_idx ++;
            variables[asub[0]].emplace_back(constraint_index, aval[0]);
            variables[asub[1]].emplace_back(constraint_index, aval[1]);
            variables[asub[2]].emplace_back(constraint_index, aval[2]);
            constraint_index ++;
        }
    }

    /*   joint points  */
    //ROS_WARN(" joint position");
    {
        int sub_shift = 0;
        double val0, val1;
        for(int k = 0; k < (segment_num - 1); k ++ )
        {   
            double scale_k = corridor[k].t;
            double scale_n = corridor[k+1].t;
            // position :
            val0 = scale_k;
            val1 = scale_n;
            for(int i = 0; i < 3; i++)
            {  // loop for x, y, z
                int nzi = 2;
                int asub[nzi];
                double aval[nzi];

                // This segment's last control point
                aval[0] = 1.0 * val0;
                asub[0] = sub_shift + (i+1) * s1d1CtrlP_num - 1;

                // Next segment's first control point
                aval[1] = -1.0 * val1;
                asub[1] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num;

                variables[asub[0]].emplace_back(constraint_index, aval[0]);
                variables[asub[1]].emplace_back(constraint_index, aval[1]);
                //variables[asub[2]].pushback(constraint_index, aval[2]);
                constraint_index ++;
                // r = MSK_putarow(task, row_idx, nzi, asub, aval);    
                // row_idx ++;
            }
            
            for(int i = 0; i < 3; i++)
            {  
                int nzi = 4;
                int asub[nzi];
                double aval[nzi];
                
                // This segment's last velocity control point
                aval[0] = -1.0;
                aval[1] =  1.0;
                asub[0] = sub_shift + (i+1) * s1d1CtrlP_num - 2;    
                asub[1] = sub_shift + (i+1) * s1d1CtrlP_num - 1;   
                // Next segment's first velocity control point
                aval[2] =  1.0;
                aval[3] = -1.0;

                asub[2] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num;    
                asub[3] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num + 1;

                variables[asub[0]].emplace_back(constraint_index, aval[0]);
                variables[asub[1]].emplace_back(constraint_index, aval[1]);
                variables[asub[2]].emplace_back(constraint_index, aval[2]);
                variables[asub[3]].emplace_back(constraint_index, aval[3]);
                constraint_index ++;
            }
            // acceleration :
            val0 = 1.0 / scale_k;
            val1 = 1.0 / scale_n;
            for(int i = 0; i < 3; i++)
            {  
                int nzi = 6;
                int asub[nzi];
                double aval[nzi];
                
                // This segment's last velocity control point
                aval[0] =  1.0  * val0;
                aval[1] = -2.0  * val0;
                aval[2] =  1.0  * val0;
                asub[0] = sub_shift + (i+1) * s1d1CtrlP_num - 3;    
                asub[1] = sub_shift + (i+1) * s1d1CtrlP_num - 2;   
                asub[2] = sub_shift + (i+1) * s1d1CtrlP_num - 1;   
                // Next segment's first velocity control point
                aval[3] =  -1.0  * val1;
                aval[4] =   2.0  * val1;
                aval[5] =  -1.0  * val1;
                asub[3] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num;    
                asub[4] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num + 1;
                asub[5] = sub_shift + s1CtrlP_num + i * s1d1CtrlP_num + 2;

                variables[asub[0]].emplace_back(constraint_index, aval[0]);
                variables[asub[1]].emplace_back(constraint_index, aval[1]);
                variables[asub[2]].emplace_back(constraint_index, aval[2]);
                variables[asub[3]].emplace_back(constraint_index, aval[3]);
                variables[asub[4]].emplace_back(constraint_index, aval[4]);
                variables[asub[5]].emplace_back(constraint_index, aval[5]);
                constraint_index ++;
            }

            sub_shift += s1CtrlP_num;
        }
    }

    for(int k = 0; k < segment_num; k++)
    {   
        Cube cube_     = corridor[k];
        double scale_k = cube_.t;

        for(int i = 0; i < 3; i++ )
        {   
            for(int j = 0; j < n_poly; j ++ )
            {   
                // pair<MSKboundkeye, pair<double, double> > vb_x;
                // pair<double, double> vb_x;
                int asub;
                double aval;
                double lo_bound, up_bound;
                asub = k * s1CtrlP_num + i * s1d1CtrlP_num + j;
                aval = 1.0;
                if(k > 0)
                {
                    lo_bound = (cube_.box[i].first  + margin) / scale_k;
                    up_bound = (cube_.box[i].second - margin) / scale_k;
                }
                else
                {
                    lo_bound = (cube_.box[i].first)  / scale_k;
                    up_bound = (cube_.box[i].second) / scale_k;
                }

                // vb_x  = make_pair( lo_bound, up_bound ) ; // # vb_x means: varialbles boundary of unknowns x (Polynomial coeff)
                variables[asub].emplace_back(constraint_index, aval);
                lower_bounds.push_back(lo_bound);
                upper_bounds.push_back(up_bound);
                // var_bdk.push_back(vb_x);
                constraint_index++;
            }
        } 
    }

    int ind_p = 0;
    for (int i = 0; i < ctrlP_num; ++i) {
        A_indptr.push_back(ind_p);
        for (const auto& variable_nz : variables[i]) {
            // coefficient
            A_data.push_back(variable_nz.second);

            // constraint index
            A_indices.push_back(variable_nz.first);
            ++ind_p;
        }
    }
    // We indeed need this line because of
    // https://github.com/oxfordcontrol/osqp/blob/master/src/cs.c#L255
    A_indptr.push_back(ind_p);


    //ROS_WARN("[Bezier Trajectory] Start stacking the objective");
    
    int min_order_l = floor(minimize_order);
    int min_order_u = ceil (minimize_order);

    int NUMQNZ = 0;
    for(int i = 0; i < segment_num; i ++)
    {   
        int NUMQ_blk = (traj_order + 1);                       // default minimize the jerk and minimize_order = 3
        NUMQNZ      += 3 * NUMQ_blk * (NUMQ_blk + 1) / 2;
    }

    vector<c_float> P_data;
    vector<c_int> P_indices;
    vector<c_int> P_indptr;
    {    
        int sub_shift = 0;
        int idx = 0;
        for(int k = 0; k < segment_num; k ++)
        {
            double scale_k = corridor[k].t;
            for(int p = 0; p < 3; p ++ )
            {
                for( int j = 0; j < s1d1CtrlP_num; j ++ )
                {
                    P_indptr.push_back(idx);
                    for( int i = 0; i < s1d1CtrlP_num; i ++ )
                    {
                        if( j >= i )
                        {
                            P_indices.push_back(sub_shift + p * s1d1CtrlP_num + i);
                            // qsubi[idx] = sub_shift + p * s1d1CtrlP_num + i; 
                            // qsubj[idx] = sub_shift + p * s1d1CtrlP_num + j;  
                            //qval[idx]  = MQM(i, j) /(double)pow(scale_k, 3);
                            if(min_order_l == min_order_u)
                                P_data.push_back(MQM(i, j) /(double)pow(scale_k, 2 * min_order_u - 3));
                                // qval[idx]  = MQM(i, j) /(double)pow(scale_k, 2 * min_order_u - 3);
                                
                            else
                                P_data.push_back(( (minimize_order - min_order_l) / (double)pow(scale_k, 2 * min_order_u - 3)
                                              + (min_order_u - minimize_order) / (double)pow(scale_k, 2 * min_order_l - 3) ) * MQM(i, j));
                                // qval[idx] = ( (minimize_order - min_order_l) / (double)pow(scale_k, 2 * min_order_u - 3)
                                            // + (min_order_u - minimize_order) / (double)pow(scale_k, 2 * min_order_l - 3) ) * MQM(i, j);
    
                            idx ++ ;
                        }
                    }
                }
            }

            sub_shift += s1CtrlP_num;
        }
        P_indptr.push_back(idx);
    }
    vector<c_float> q;
    for (int i = 0; i < ctrlP_num; i++)
    {
        q.push_back(0.0);
    }


    // OSQPData* data = FormulateProblem();
    OSQPData* data = reinterpret_cast<OSQPData*>(c_malloc(sizeof(OSQPData)));
    // CHECK_EQ(lower_bounds.size(), upper_bounds.size());

    // size_t kernel_dim = 3 * num_of_knots_;
    // size_t num_affine_constraint = lower_bounds.size();
    size_t kernel_dim = ctrlP_num;
    size_t num_affine_constraint = lower_bounds.size();

    data->n = kernel_dim;
    data->m = num_affine_constraint;
    data->P = csc_matrix(kernel_dim, kernel_dim, NUMQNZ, CopyData(P_data),
                         CopyData(P_indices), CopyData(P_indptr));
    //P_data.size()
    data->q = CopyData(q);
    data->A = csc_matrix(num_affine_constraint, kernel_dim, A_data.size(),
                   CopyData(A_data), CopyData(A_indices), CopyData(A_indptr));
    data->l = CopyData(lower_bounds);
    data->u = CopyData(upper_bounds);

    ros::Time time_end1 = ros::Time::now();


    OSQPSettings* settings = SolverDefaultSettings();
    // settings->max_iter = max_iter;
    OSQPWorkspace* osqp_work = osqp_setup(data, settings);
    osqp_solve(osqp_work);
    auto status = osqp_work->info->status_val;

    if (status < 0 || (status != 1 && status != 2)) {
        osqp_cleanup(osqp_work);
        FreeData(data);
        c_free(settings);
        return false;
    } else if (osqp_work->solution == nullptr) {
        osqp_cleanup(osqp_work);
        FreeData(data);
        c_free(settings);
        return false;
    }

    ros::Time time_end2 = ros::Time::now();
    ROS_WARN("time consume in optimize is :");
    cout<<time_end2 - time_end1<<endl;


    VectorXd d_var(ctrlP_num);
    for(int i = 0; i < ctrlP_num; i++)
        d_var(i) = osqp_work->solution->x[i];
    
    PolyCoeff = MatrixXd::Zero(segment_num, 3 *(traj_order + 1) );

    int var_shift = 0;
    for(int i = 0; i < segment_num; i++ )
    {
        for(int j = 0; j < 3 * n_poly; j++)
            PolyCoeff(i , j) = d_var(j + var_shift);

        var_shift += 3 * n_poly;
    }   

    // Cleanup
    osqp_cleanup(osqp_work);
    FreeData(data);
    c_free(settings);
    return true;
}

OSQPSettings* TrajectoryGenerator::SolverDefaultSettings() {
    // Define Solver default settings
    OSQPSettings* settings =
        reinterpret_cast<OSQPSettings*>(c_malloc(sizeof(OSQPSettings)));
    osqp_set_default_settings(settings);
    settings->polish = true;
    settings->verbose = false;
    settings->scaled_termination = true;
    return settings;
}
void TrajectoryGenerator::FreeData(OSQPData* data) {
    delete[] data->q;
    delete[] data->l;
    delete[] data->u;

    delete[] data->P->i;
    delete[] data->P->p;
    delete[] data->P->x;

    delete[] data->A->i;
    delete[] data->A->p;
    delete[] data->A->x;
}