#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  if(estimations.size() != ground_truth.size() || estimations.size() == 0)
      return rmse;

  //accumulate squared residuals
  VectorXd error;
  for(int i=0; i < estimations.size(); ++i){
      error = estimations[i]-ground_truth[i];
      error = error.array()*error.array();
      rmse = rmse + error;
  }

  //calculate the mean
  rmse = rmse / estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);
	float pxpy = px*px+py*py;
	//check division by zero
	if(pxpy < 0.0001)
	{
    cout<<"CalculateJacobian() - Error - Division by Zero"<<endl;
    pxpy = 0.0001;
    //return Hj; //is that a good idea? probably not.
  }
	float sqpxpy = sqrt(pxpy);
	float pxpy32 = pxpy*sqpxpy;
	float vpvp = vx*py-vy*px;

  Hj(0,0) = px/sqpxpy;
  Hj(0,1) = py/sqpxpy;
  Hj(0,2) = 0;
  Hj(0,3) = 0;
  Hj(1,0) = -py/pxpy;
  Hj(1,1) = px/pxpy;
  Hj(1,2) = 0;
  Hj(1,3) = 0;
  Hj(2,0) = py*vpvp/pxpy32;
  Hj(2,1) = px*vpvp/pxpy32;
  Hj(2,2) = px/sqpxpy;
  Hj(2,3) = py/sqpxpy;

	return Hj;
}
