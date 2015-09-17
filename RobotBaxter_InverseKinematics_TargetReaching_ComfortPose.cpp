/*
 *  RajeevRatan_S1471378.cpp
 *
 *  Created on: 25 Feb 2014
 *      Author: Rajeev Ratan
 *      Demonstration of Baxter Arm Control - Inverse Kinematics - Target Rearching using direct Jacobian, then with a comfort pose
 */

#include "BaxterTools.h"
#include "unistd.h"
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

Eigen::MatrixXd J_PSEUDO (char A, Eigen::MatrixXd J, Eigen::MatrixXd Winv, Eigen::MatrixXd Cinv);
Eigen::VectorXd Drift (int comf,char A,Eigen::VectorXd q);
Eigen::VectorXd Delta_Y (char A,int target_num, Eigen::VectorXd y,Eigen::VectorXd target, Eigen::VectorXd ystart, int itr);
Eigen::VectorXd Update_Q (char A, Eigen::VectorXd q, Eigen::VectorXd delta_y, Eigen::VectorXd h, Eigen::MatrixXd J, Eigen::MatrixXd J_Hash, Eigen::MatrixXd I, int target_num, int no_null,char smooth);
float Error (char A, int target_num, Eigen::VectorXd y,Eigen::VectorXd target);
float WCost (char A,int start,int target_num,Eigen::VectorXd q,Eigen::MatrixXd W);
void Write_Text_Q1 (vector<double> question1);
Eigen::VectorXd Get_Start(int start);


int main(int argc,char* argv[])
{
  std::cout << "\033[2J\033[1;1H"; //Clear Terminal

  std::cout << "RLSC Assignment\n";
  std::cout << "Rajeev Ratan S1471378\n";
  std::cout << "\nProceeding through Questions 1,2,3 and 5\n";
  std::cout << "\n\nPlease wait for Baxter to get his act together....should take about 10 seconds.\n";

  // Create the robot interface object
  BaxterTools bax;

  // Connect to the simulator
  if(argc==2)
    {
      bax.Connect(argv[1]);
    }
  else
    {
      bax.Connect("localhost");
    }
  // Start the simulation
  bax.StartSimulation();

  double t=0;

  Eigen::VectorXd q=Eigen::VectorXd::Zero(18); // Joint angles
  Eigen::VectorXd y; // End-effector position
  Eigen::VectorXd target; // Target positions
  Eigen::MatrixXd J; // Jacobian matrix

  //////////////////////////////////////////////////////////////////////
  Eigen::MatrixXd Winv = Eigen::MatrixXd::Zero(7,7); // Weighting matrix
  for(int i=0;i<7;i++) Winv(i,i) = ((double)i)/6.0+0.1;

  Eigen::MatrixXd C = Eigen::MatrixXd::Identity(3,3)*1e3; // Regularisation
  Eigen::MatrixXd Cinv = Eigen::MatrixXd::Identity(3,3)*1e-3;

  Eigen::VectorXd qstart1(18); // Starting pose 1
  Eigen::VectorXd qstart2(18); // Starting pose 2
  Eigen::VectorXd qstart3(18); // Starting pose 3
  qstart1 << M_PI/4.0,0,0,M_PI/2.0,0,0,0,0,0,      -M_PI/4.0,0,0,M_PI/2.0,0,0,0,0,0;
  qstart2 << -M_PI/4.0,0,0,M_PI/2.0,M_PI/2.0,M_PI/2.0,0,0,0,   M_PI/4.0,0,0,M_PI/2.0,-M_PI/2.0,M_PI/2.0,0,0,0;
  qstart3 << M_PI/4.0,-M_PI/2.0,0,M_PI/2.0,-M_PI/4.0,-M_PI/4.0,0,0,0,      -M_PI/4.0,-M_PI/2.0,0,M_PI/2.0,M_PI/4.0,-M_PI/4.0,0,0,0;

  Eigen::VectorXd q_comf1(18); // Comfortable pose 1
  Eigen::VectorXd q_comf2(18); // Comfortable pose 2
  q_comf1 << 0,-0.5,0,0.5,0,1.5,0,0,0,  0,-0.5,0,0.5,0,1.5,0,0,0 ;
  q_comf2 << -20.0/180.0*M_PI, 40.0/180.0*M_PI, 70.0/180.0*M_PI, 90.0/180.0*M_PI, 0,0.5,0,0,0, 20.0/180.0*M_PI, 40.0/180.0*M_PI, -70.0/180.0*M_PI, 90.0/180.0*M_PI, 0,0.5,0,0,0;

  //////////////////////////////////////////////////////////////////////

  // More Variables & Constants
  Eigen::MatrixXd I = Eigen::Matrix<double, 7, 7>::Identity();
  Eigen::VectorXd qold=Eigen::VectorXd::Zero(18);
  Eigen::MatrixXd Jinv; // Jacobian Pseudo Inverse matrix
  Eigen::VectorXd ytarget;
  Eigen::VectorXd delta_y;
  Eigen::VectorXd new_q;
  Eigen::VectorXd h;
  Eigen::VectorXd ystart;
  int k = 0;
  bax.GetTargets(target);
  float error = 100;
  Eigen::MatrixXd W = Winv.inverse();
  double costs; 
  int target_num;  // Set Target
  int itr;
  char key=0;
  int c = 0;

  //Set defaults
  char A = 'R';             // Set Arm
  int comf = 1;             // Set Comfort Matrix
  int  start = 1;           // Set Starting Position
  int  p = 0;               // Target Count
  int n = 0;
  Eigen::VectorXd store;
  vector<double> question1;
  int num = 0;
  int s = 0;
  int no_null = 0;
  float E = 0.02;
  char smooth = 0.05;

  std::cout << "Moving one arm at a time" << "\n\n";
 
  while(key!='q')
    {
      while (p<8){    //Iterate through targets 1 to 8;
	error = 100;
	no_null = 0;
	q = qstart1;
	bax.SetJointAngles(q);
	sleep(1);

	k = 0;
	ystart=bax.GetIK(q);

	int itr = 1; //test

	while (error > E){

	  // Get pressed key
	  key=bax.GetKey();

	  // My Code ////////////////////
	  k++;
	  y=bax.GetIK(q);      // Get end-effector position
	  J=bax.GetJ(q);       // Get Jacobian of the end effector+
	  target_num = p;      // Set Target
	  delta_y = Delta_Y(A,target_num,y,target,ystart,itr);
	  h = Drift(comf,A,q);
	  Eigen::MatrixXd J_Hash = J_PSEUDO(A,J,Winv,Cinv);
	  q = Update_Q(A,q,delta_y,h,J,J_Hash,I,target_num,no_null,smooth);

	  bax.SetJointAngles(q);
	  bax.AdvanceSimulation();	 
	  y=bax.GetIK(q);
	  error = Error(A,target_num,y,target);
	}
	
	costs = WCost(A,start,target_num,q,W);
	std::cout << "Using Arm - " << A << " and Comfort Position " << comf-1 << " for Target # " << p << "\n";
	std::cout << k << " Iterations to Final Q, " <<"Cost = " << costs << "\n";

	question1.push_back(costs);
	break;
      }
      p++; 
          
      if (p==8){   

	std::cout << "Cycle = " << c+1 << " Completed\n";

	if( c == 0){
	  p = 0;                          // restart target count from 0
	  if(A == 'L'){A='R';}            // switch arms
	  else{A ='L';}
	  std::cout << "\nChanging Arms" << "\n\n";
	}

	if (c == 1){                    // Switch comfort poses after reaching target with both arms
	  if(A == 'R'){A='L';}          // Switch arms
	  else{A ='R';}
	  std::cout << "\nChanging Arms" << "\n\n";
	  comf = 2;                    //Change comfort position to 2
	  p = 0;                       //restart target count from 0
	  std::cout << "Changing to comfort pose " << comf-1 << "\n\n";
	}

	if (c == 2){ 
	  if(A == 'L'){A='R';}        
	  else{A ='L';}
	  p = 0; //restart target count from 0
	  std::cout << "\nChanging Arms" << "\n\n";
	  comf = 2;                    //Change comfort position to 2 
	}

	if (c == 3){
	  if(A = 'R'){A='L';}// switch arms
	  else{A ='R';}
	  p = 0; //restart target count from 0
	  std::cout << "\nChanging Arms" << "\n\n";
	}

	c++;

	if(c == 4){
	  Write_Text_Q1 (question1);
	  std::cout << "Question 1 Completed" << "\n\n";
	  num = 2;
	}
      }	  
      
      if (num == 2){
 

	//Set Parameters
	A = 'R';             // Set Arm
	comf = 1;            // Set Comfort Matrix
	start = 1;           // Set Starting Position
	p = 0;               // Target Count
	k = 0;
	error = 100;
	string Comfort;
	no_null = 1;

	std::cout << "Using Right Arm for reaching Target 0 with Comfort Position " << comf-1 << "\n\n";	

	while (start < 4){

	  q=Get_Start(start);
	  bax.SetJointAngles(q);
	  bax.AdvanceSimulation();

	  error = 100;
	  k++;
	  int iterations = 0;

	  while (error > E){
	
	    y=bax.GetIK(q);      // Get end-effector position
	    J=bax.GetJ(q);       // Get Jacobian of the end effector+
	    target_num = p;      // Set Target
	    delta_y = Delta_Y(A,target_num,y,target,ystart,itr);
	    h = Drift(comf,A,q);
	    if(k > 3){no_null = 0;}
	    Eigen::MatrixXd J_Hash = J_PSEUDO(A,J,Winv,Cinv);
	    q = Update_Q(A,q,delta_y,h,J,J_Hash,I,target_num,no_null,smooth);
	    bax.SetJointAngles(q);
	    bax.AdvanceSimulation();
	    y=bax.GetIK(q);
	    error = Error(A,target_num,y,target);
	    iterations++;
	  }

	  if(no_null == 1){ Comfort = "Off";}
	  else{Comfort = "On";}
	  std::cout <<"Right Arm, Target 0, " << "Starting Position Q" << start << ", Comfort Pose 0 - " << Comfort << "\n";
	  std::cout << "Took " << iterations << " Iterations to get to Target\n";
	  // sleep(1);
	  start++;
	  if(start == 4 && k == 3){start = 1;}
	}


	//Set Parameters
	A = 'R';             // Set Arm
	comf = 1;            // Set Comfort Matrix
	start = 1;           // Set Starting Position
	p = 0;               // Target Count
	k = 0;
	error = 100;
	no_null = 0;

	std::cout << "Reaching all 8 Targets with Right Arm for each of the 3 starting pose using Comfort Pose " << comf-1 << "\n\n";	

	ofstream myfile;
	myfile.open ("matrix.txt");
	myfile << "Question 3 - Matrix\n \n";
	Eigen::MatrixXd Question3;

	while (start < 4){

	  q=Get_Start(start);
	  bax.SetJointAngles(q);
	  bax.AdvanceSimulation();

	  error = 100;
	  k++;

	  while (error > E){
	
	    y=bax.GetIK(q);      // Get end-effector position
	    J=bax.GetJ(q);       // Get Jacobian of the end effector+
	    target_num = p;      // Set Target
	    delta_y = Delta_Y(A,target_num,y,target,ystart,itr);
	    h = Drift(comf,A,q);

	    Eigen::MatrixXd J_Hash = J_PSEUDO(A,J,Winv,Cinv);
	    q = Update_Q(A,q,delta_y,h,J,J_Hash,I,target_num,no_null,smooth);
	    bax.SetJointAngles(q);
	    bax.AdvanceSimulation();
	    y=bax.GetIK(q);
	    error = Error(A,target_num,y,target);
	  }

	  myfile << q.transpose() << "\n";

	  if(no_null == 1){ Comfort = "Off";}
	  else{Comfort = "On";}
	  std::cout <<"Target " << p << ", Starting Position Q_start" << start << ", Comfort Pose 0 - " << Comfort << "\n";

	  start++;
	  if(start == 4){start = 1;p++;}
	  if(p == 8){
	    myfile.close();
	    std::cout << "\n24 x 7 Matrix created in file matrix.txt located in the root directory\n";
	    break;}
	}
      


	//Set Parameters
	A = 'B';             // Set Arm
	comf = 1;            // Set Comfort Matrix
	start = 1;           // Set Starting Position
	p = 7;               // Target 
	k = 0;
	error = 100;
	float error2 = 100;
	float error_ave = 100;
	no_null = 0;
	float E = 0.01;
	smooth = 1; //Smooth Motion turned on
	Eigen::MatrixXd J_Hash;
	std::cout << "Moving both arms from Starting Position 1 to Target 7 using comfort pose 0" << "\n";

	q=Get_Start(start);
	int exit = 0;

	bax.SetJointAngles(q);
	bax.AdvanceSimulation();

	while (error && error2 > E){

	  y=bax.GetIK(q);      // Get end-effector position
	  J=bax.GetJ(q);       // Get Jacobian of the end effector+
	  target_num = p;      // Set Target

	  A ='L';
	  delta_y = Delta_Y(A,target_num,y,target,ystart,itr);
	  h = Drift(comf,A,q);
	  J_Hash = J_PSEUDO(A,J,Winv,Cinv);
	  q = Update_Q(A,q,delta_y,h,J,J_Hash,I,target_num,no_null,smooth);
	  error = Error(A,target_num,y,target);
	    
	  A ='R';
	  delta_y = Delta_Y(A,target_num,y,target,ystart,itr);
	  h = Drift(comf,A,q);
	  J_Hash = J_PSEUDO(A,J,Winv,Cinv);
	  q = Update_Q(A,q,delta_y,h,J,J_Hash,I,target_num,no_null,smooth);
	  bax.SetJointAngles(q);
	  bax.AdvanceSimulation();
	  y=bax.GetIK(q);
	  error2 = Error(A,target_num,y,target);
	}

	std::cout <<"Final Error in Left Arm = " << error << "\n";
	std::cout <<"Final Error in Right Arm = " << error2 << "\n";

	//Exit Program
	std::cout << "\nDemonstration Completed. Exiting Program!\n";
	break;
      }

    }
  
  bax.StopSimulation();
  return(0);
}



//**************** FUNCTIONS ************

//Returns J#
Eigen::MatrixXd J_PSEUDO (char A, Eigen::MatrixXd J, Eigen::MatrixXd Winv, Eigen::MatrixXd Cinv)
{
  Eigen::MatrixXd J_Pos;
  if (A == 'R'){
    J_Pos = J.block(0,0,3,7);
  }
  else if(A == 'L'){
    J_Pos = J.block(6,7,3,7);
  }
  else{std::cout << "No correct arm selected!\n";
  }
  //J Transpose
  Eigen::MatrixXd J_Transpose = J_Pos.transpose(); 
  //Calculating the inside elements of the J# 
  Eigen::MatrixXd JS2 = (((J_Pos * Winv) * J_Transpose) + Cinv);
  //Inverse of the inside brackets
  Eigen::MatrixXd JSinv = JS2.inverse();
  //(Winv * J_Tranpose) * Inside Bracket
  Eigen::MatrixXd J_PSEUDO = (Winv * J_Transpose) * JSinv;
  return J_PSEUDO;
}


//Returns h
Eigen::VectorXd Drift (int comf,char A,Eigen::VectorXd q)
{
  Eigen::VectorXd q_comf1(18); // Comfortable pose 1
  Eigen::VectorXd q_comf2(18); // Comfortable pose 2
  q_comf1 << 0,-0.5,0,0.5,0,1.5,0,0,0,  0,-0.5,0,0.5,0,1.5,0,0,0 ;
  q_comf2 << -20.0/180.0*M_PI, 40.0/180.0*M_PI, 70.0/180.0*M_PI, 90.0/180.0*M_PI, 0,0.5,0,0,0, 20.0/180.0*M_PI, 40.0/180.0*M_PI, -70.0/180.0*M_PI, 90.0/180.0*M_PI, 0,0.5,0,0,0;

  int arm;
  if (A == 'R'){arm = 0;}
  else if(A == 'L'){arm = 9;}
  else{std::cout << "No correct arm selected!\n";
  }
 
  //Select Comfort Matrix
  if (comf == 1){
    return (q_comf1.segment(arm,7) - q.segment(arm,7));    
  }
  if (comf == 2){
    return (q_comf2.segment(arm,7) - q.segment(arm,7));    
  }  
}


//Returns (Y_Target - Current_Y_Position)
Eigen::VectorXd Delta_Y (char A,int target_num, Eigen::VectorXd y,Eigen::VectorXd target, Eigen::VectorXd ystart, int itr)
{
  int arm;
  Eigen::VectorXd dy;
  int t = target_num * 3;

  if (A == 'R'){
    arm = 0;
    return (target.segment(t,3) - y.segment(arm,3));
 
  }
  else if(A == 'L'){
    arm = 6;   
    return (target.segment(t,3) - y.segment(arm,3));
  }
  else{std::cout << "No correct arm selected!\n";
  }
}


//Updates Q with dq
Eigen::VectorXd Update_Q (char A, Eigen::VectorXd q, Eigen::VectorXd delta_y, Eigen::VectorXd h, Eigen::MatrixXd J, Eigen::MatrixXd J_Hash, Eigen::MatrixXd I, int target_num, int no_null,char smooth)
{
  // int arm;
  int armq;
  Eigen::MatrixXd J_Pos;
  Eigen::VectorXd Null;
  Eigen::VectorXd DQ;
  Eigen::VectorXd tq;
  float scaler = 0.3;

  if (smooth == 1){scaler = 0.05;}

  if (A == 'R'){
    //  arm = 0;
    armq = 0;
    J_Pos = J.block(0,0,3,7);
  }
  else if(A == 'L'){
    //  arm = 6;
    armq = 9;
    J_Pos = J.block(6,7,3,7);
  }
  else{std::cout << "No correct arm selected!\n";}

  Null =  ((I - (J_Hash*J_Pos)) * h);

  if(no_null == 1){
    DQ =  (J_Hash * delta_y);
    tq = q.segment(armq,7) + (DQ*scaler);
  }
  else{
    DQ =  (J_Hash * delta_y) + Null;
    tq = q.segment(armq,7) + (DQ*scaler);
  }
  
  for(int i=0;i<7;i++){ q(i+armq)=tq(i); }
  return q;
}


//Calculates the Position Error of the End Effector from the Target
float Error (char A, int target_num, Eigen::VectorXd y,Eigen::VectorXd target)
{
  int arm;
  int t = target_num * 3;
  Eigen::VectorXd Tar;
  Eigen::VectorXd EndEffectorPos;
  Eigen::VectorXd Dist;
  if (A == 'R'){ arm = 0;}
  else if(A == 'L'){arm = 6;}
  else{std::cout << "No correct arm selected!\n";}

  EndEffectorPos = y.segment(arm,3);
  Tar = target.segment(t,3);
  Dist = Tar - EndEffectorPos;
  return Dist.norm();
}


//Calculates the Cost (Q_Final - Q_Start)
float WCost (char A,int start,int target_num,Eigen::VectorXd q,Eigen::MatrixXd W)
{
  int arm;
  int error; 
  int t = target_num * 3;
  Eigen::VectorXd q_delta;
  Eigen::VectorXd qstart(18);  // Starting pose 1
  Eigen::VectorXd qstart1(18); // Starting pose 1
  Eigen::VectorXd qstart2(18); // Starting pose 2
  Eigen::VectorXd qstart3(18); // Starting pose 3
  qstart1 << M_PI/4.0,0,0,M_PI/2.0,0,0,0,0,0,      -M_PI/4.0,0,0,M_PI/2.0,0,0,0,0,0;
  qstart2 << -M_PI/4.0,0,0,M_PI/2.0,M_PI/2.0,M_PI/2.0,0,0,0,   M_PI/4.0,0,0,M_PI/2.0,-M_PI/2.0,M_PI/2.0,0,0,0;
  qstart3 << M_PI/4.0,-M_PI/2.0,0,M_PI/2.0,-M_PI/4.0,-M_PI/4.0,0,0,0,      -M_PI/4.0,-M_PI/2.0,0,M_PI/2.0,M_PI/4.0,-M_PI/4.0,0,0,0;

  if(start == 1){qstart = qstart1;}
  if(start == 2){qstart = qstart2;}
  if(start == 3){qstart = qstart3;}

  if (A == 'R'){ arm = 0;}
  else if(A == 'L'){arm = 9;}
  else{std::cout << "No correct arm selected!\n";}

  q_delta = qstart.segment(arm,7) - q.segment(arm,7); 
  return (q_delta.transpose() * W) * q_delta;
}


//Creates output text file with table for Question 1
void Write_Text_Q1 (vector<double> question1)
{
  ofstream myfile;
  myfile.open ("output.txt");
  std::cout << "Writing to File\n";
  myfile << "Question 1 - Cost Table Output\n \n";
  myfile << "\tComfort Pose 1\t\t\tComfort Pose 2 \n";
  myfile << "Target\tRight\tLeft\tBest Arm\tRight\tLeft\tBest Arm\n";

  for(int l=0;l<8;l++){ 
    myfile << l << "\t" << question1[l] << "\t" <<question1[l+8];
    if(question1[l]<question1[l+8]){
      myfile << "\tRight"; }
    else{
      myfile << "\tLeft";
    }
    myfile << "\t\t" << question1[l+16] << "\t" << question1[l+24];
    if(question1[l+16]<question1[l+24]){
      myfile << "\tRight\n"; }
    else{
      myfile << "\tLeft\n";
    }     
  }
  myfile.close();
}


//Gets Starting Position
Eigen::VectorXd Get_Start(int start)
{
  Eigen::VectorXd qstart(18);  // Starting pose 1
  Eigen::VectorXd qstart1(18); // Starting pose 1
  Eigen::VectorXd qstart2(18); // Starting pose 2
  Eigen::VectorXd qstart3(18); // Starting pose 3
  qstart1 << M_PI/4.0,0,0,M_PI/2.0,0,0,0,0,0,      -M_PI/4.0,0,0,M_PI/2.0,0,0,0,0,0;
  qstart2 << -M_PI/4.0,0,0,M_PI/2.0,M_PI/2.0,M_PI/2.0,0,0,0,   M_PI/4.0,0,0,M_PI/2.0,-M_PI/2.0,M_PI/2.0,0,0,0;
  qstart3 << M_PI/4.0,-M_PI/2.0,0,M_PI/2.0,-M_PI/4.0,-M_PI/4.0,0,0,0,      -M_PI/4.0,-M_PI/2.0,0,M_PI/2.0,M_PI/4.0,-M_PI/4.0,0,0,0;

  if(start == 1){return qstart1;}
  if(start == 2){return qstart2;}
  if(start == 3){return qstart3;}
}

