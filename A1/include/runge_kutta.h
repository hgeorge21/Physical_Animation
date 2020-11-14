//Input:
//  q - generalized coordiantes for the mass-spring system
//  qdot - generalized velocity for the mass spring system
//  dt - the time step in seconds
//  mass - the mass
//  force(q, qdot) - a function that computes the force acting on the mass as a function. This takes q and qdot as parameters.
//Output:
//  q - set q to the updated generalized coordinate using Runge-Kutta time integration
//  qdot - set qdot to the updated generalized velocity using Runge-Kutta time integration

template<typename FORCE> 
inline void runge_kutta(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt, double mass,  FORCE &force) {
    Eigen::VectorXd k1(q.size());
    Eigen::VectorXd k2(q.size());
    Eigen::VectorXd k3(q.size());
    Eigen::VectorXd k4(q.size());

    Eigen::VectorXd temp;
    force(k1, q, qdot);

    temp = q+dt*(k1/mass);
    force(k2, temp, qdot);

    temp = q+dt*(k2/mass);
    force(k3, temp, qdot);

    temp = q+dt*(k3/mass);
    force(k4, temp, qdot);


    qdot = qdot + dt/6*(k1+2*k2+2*k3+k4)/mass;
    q = q + dt*qdot;

}