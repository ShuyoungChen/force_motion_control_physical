import general_robotics_toolbox as rox
import numpy as np
import threading
import quadprog
from scipy.linalg import norm

class QuadProg(object):
	def __init__(self, robot):	
		self._robot = robot  
		self._robot.joint_vel_limit = np.deg2rad(np.array([110, 90, 90, 190, 140, 235]))
		self._lock = threading.Lock()
		
		# inequality constraints for quadprog
		self._h = np.zeros((12, 1))
		self._sigma = np.zeros((12, 1))
		self._dhdq = np.vstack((np.hstack((np.eye(6), np.zeros((6, 1)), np.zeros((6, 1)))), np.hstack((-np.eye(6), np.zeros((6, 1)), np.zeros((6, 1))))))

		# parameters for solving QP
		self._c = 0.5
		self._eta = 0.1
		self._E = 0.0005
		self._Ke = 1

		# optimization params to formalize object function
		self._er = 0.01
		self._ep = 0.01
		self._epsilon = 0.05
		
		# number of joints
		self._n = len(self._robot.joint_type)
		
		# for fixed orientation equality constraint
		self._eq_orien = False
		
		# desired orientation of eef
		self._R_des = np.eye(3)
		
	def compute_joint_vel_cmd_qp(self, joint_position, spatial_velocity_command):		
		J_eef = rox.robotjacobian(self._robot, joint_position)

		# desired rotational velocity
		vr = spatial_velocity_command[0:3]
		vr = vr.reshape(3, 1)
		
		# desired linear velocity
		vp = spatial_velocity_command[3:None]
		vp = vp.reshape(3, 1)	
		
		Q = self.getqp_H(J_eef, vr, vp)         
	
		# make sure Q is symmetric
		Q = 0.5*(Q + Q.T)
		
		f = self.getqp_f()
		f = f.reshape((8, ))

		# change the velocity scale
		LB = np.vstack((-self._robot.joint_vel_limit.reshape(6, 1), 0, 0))
		UB = np.vstack((self._robot.joint_vel_limit.reshape(6, 1), 1, 1))
	
		# inequality constrains A and b
		self._h[0:6] = joint_position.reshape(6, 1) - self._robot.joint_lower_limit.reshape(6, 1)
		self._h[6:12] = self._robot.joint_upper_limit.reshape(6, 1) - joint_position.reshape(6, 1)
		
		self._sigma[0:12] = self.inequality_bound(self._h[0:12])
		
		A = np.vstack((self._dhdq, np.eye(8), -np.eye(8)))
		b = np.vstack((self._sigma, LB, -UB))
		b = b.reshape((28, ))
		
		# solve the quadprog problem
		# scale the matrix to avoid numerical errors of solver
		sc = norm(Q,'fro')
		dq_sln = quadprog.solve_qp(Q/sc, -f/sc, A.T, b)[0]
			
		if len(dq_sln) < self._n:
			qdot = np.zeros((self._n,1))
			V_scaled = 0
			print 'No Solution'
		else:
			qdot = dq_sln[0: self._n]
			V_scaled = dq_sln[-1]*vp
	
		V_linear = np.dot(J_eef[3:6,:], qdot)
		V_rot = np.dot(J_eef[0:3,:], qdot)
		
		qdot = qdot.reshape((6, ))
		
		#print 'desired angular velocity'
		#print vr
		#print 'actual angular velocity'
		#print V_rot
		#print 'desired linear velocity'
		#print vp
		#print 'actual linear velocity'
		#print V_linear
	
		return qdot
			
	# for inequality constraint        
	def inequality_bound(self, h):
		sigma = np.zeros((h.shape))
		h2 = h - self._eta
		sigma[np.array(h2 >= self._epsilon)] = -np.tan(self._c*np.pi/2)
		sigma[np.array(h2 >= 0) & np.array(h2 < self._epsilon)] = -np.tan(self._c*np.pi/2/self._epsilon*h2[np.array(h2 >= 0) & np.array(h2 < self._epsilon)])
		sigma[np.array(h >= 0) & np.array(h2 < 0)] = -self._E*h2[np.array(h >= 0) & np.array(h2 < 0)]/self._eta
		sigma[np.array(h < 0)] = self._E

		return sigma
			  
	# get the matrix f for solving QP       
	def getqp_f(self):
		f = -2*np.array([0, 0, 0, 0, 0, 0, self._er, self._ep]).reshape(8, 1)

		return f
	 
	# get the matrix H for solving QP        
	def getqp_H(self, J, vr, vp):
		H1 = np.dot(np.hstack((J,np.zeros((6,2)))).T,np.hstack((J,np.zeros((6,2)))))

		tmp = np.vstack((np.hstack((np.hstack((np.zeros((3, self._n)),vr)),np.zeros((3,1)))),np.hstack((np.hstack((np.zeros((3,self._n)),np.zeros((3,1)))),vp)))) 
		H2 = np.dot(tmp.T,tmp)

		H3 = -2*np.dot(np.hstack((J,np.zeros((6,2)))).T, tmp)
		H3 = (H3+H3.T)/2;

		tmp2 = np.vstack((np.array([0,0,0,0,0,0,np.sqrt(self._er),0]),np.array([0,0,0,0,0,0,0,np.sqrt(self._ep)])))
		H4 = np.dot(tmp2.T, tmp2)

		H = 2*(H1+H2+H3+H4)

		return H