import numpy as np
import rpi_abb_irc5
import time
import general_robotics_toolbox as rox
import QuadProg as qp
import sys
import rpi_ati_net_ft
import matplotlib.pyplot as plt

def abb_irb6640_180_255_robot():
    """Return a Robot instance for the ABB IRB6640 180-255 robot"""
    
    x = np.array([1,0,0])
    y = np.array([0,1,0])
    z = np.array([0,0,1])
    a = np.array([0,0,0])
    
    H = np.array([z,y,y,x,y,x]).T
    P = np.array([0.78*z, 0.32*x, 1.075*z, 0.2*z, 1.1425*x, 0.2*x, a]).T
    joint_type = [0,0,0,0,0,0]
    joint_min = np.deg2rad(np.array([-170, -65, -180, -300, -120, -360]))
    joint_max = np.deg2rad(np.array([170, 85, 70, 300, 120, 360]))
    
    p_tool = np.array([0,0,0])
    R_tool = rox.rot([0,1,0], np.pi/2.0)
    
    return rox.Robot(H, P, joint_type, joint_min, joint_max, R_tool=R_tool, p_tool=p_tool) 
	
def main():    
    # initialize EGM interface instance
    egm = rpi_abb_irc5.EGM()
    
    # initialize a robot instance
    abb_robot = abb_irb6640_180_255_robot()

    # desired force
    #Fd = 1000
		
    # desired velocity in y
    #vdy = -0.5
	
    # feedback gain
    Kp = 0.0002
    Kd = 0.0008
    Ki = 0.0004
	
    # time step
    delta_t = 0.004
	
    # initial configuration in degree
    init = [-91.85,2.53,38.20,0.00,49.27,-1.85]

    n = 2000

    # quadprog to solve for joint velocity
    quadprog = qp.QuadProg(abb_robot)
	
    # force sensor interface
    try:
        if (len(sys.argv) < 2):
            raise Exception('IP address of ATI Net F/T sensor required')
        host=sys.argv[1]
        netft=rpi_ati_net_ft.NET_FT(host)
        netft.set_tare_from_ft()
       
        netft.start_streaming()
    except KeyboardInterrupt:
        pass

    ####### trapezoidal desired force in z #######
    tt = np.linspace(0, 4*np.pi, n)
    desired_f = np.zeros((1, n))
    vdy = np.zeros((1, n))

    # form trap force and trap motion
    for i in range(n):
        if tt[i] >= 0 and tt[i] < np.pi:
            desired_f[0, i] = 50+302*tt[i]
            vdy[0, i] = -0.2*tt[i]
        elif tt[i] >= np.pi and tt[i] < 3*np.pi:
            desired_f[0, i] = 50+302*np.pi
            vdy[0, i] = -0.2*np.pi
        else:
            desired_f[0, i] = 50+302*np.pi-302*(tt[i]-3*np.pi)
            vdy[0, i] = -0.2*np.pi+ 0.2*(tt[i]-3*np.pi)
    #plt.plot(vdy[0, :])
    #plt.show()

    ######## change here ########
    #pos = 0
    #acce = 0
    #v_l_pre = 0
    #########
	
    # output force
    force_out = np.zeros((6, n))
	
    # pos of eef
    eef_pos = np.zeros((3, n))	
    eef_orien = np.zeros((4, n))

    # timestamp
    tim = np.zeros((1, n))

    # referece height of coupon that achieves desired force
    z_ref = 0.89226
    x_ref = 0
    y_ref = -1.35626
	
    # determine if robot has reached the initial configuration init
    tag = True
    while tag:
        res, state = egm.receive_from_robot(.1)
        if res: 
            #print sum(np.rad2deg(state.joint_angles))- sum(init)
            if np.fabs(sum(np.rad2deg(state.joint_angles)) - sum(init)) < 1e-3:
                tag = False
	
    time.sleep(1)	
    print '--------start force control--------'
	
    ### drain the force sensor buffer ###	
    for i in range(1000):	
        flag, ft = netft.read_ft_streaming(.1)
		
    ### drain the EGM buffer ###	
    for i in range(1000):	
        res, state = egm.receive_from_robot(.1)
		
    flag, ft = netft.read_ft_streaming(.1)
    F0 = ft[5]
    print F0  
    time.sleep(2) 

    cnt = 0
    while cnt < n:
        # receive EGM feedback
        res, state = egm.receive_from_robot(.1)

        if not res:
            continue
       
	    # forward kinematics
        pose = rox.fwdkin(abb_robot, state.joint_angles)
        R = pose.R
        flag, ft = netft.read_ft_streaming(.1)

        # map torque/force from sensor frame to base frame
        T_b = np.matmul(R, ft[0:3])
        F_b = np.matmul(R, ft[3:None])
        F = F_b[2] # first three torques and then three forces

        Fd0 = 50
        # do not start motion in y until robot barely touches coupon (50 N)
        if F < Fd0-0.1 and cnt == 0:
            z = pose.p[2]
            # account for the robot base and length of tool
            z = z + 0.026-0.18+0.00353
            # will shake if gain too large, here use 0.0002
            v_z = Kp*10*(F-Fd0)
            v_l = np.array([0, 0, v_z])
        else: 
            # deadzone for Fx
            if abs(F_b[0]) < 30:
                F_x = 0

            # deadzone for Fy
            if abs(F_b[1]) < 30:
                F_y = 0

            v_x = Kp/2*(F_x-0)
            v_y = vdy[0, cnt] + Kp/2*(F_y-0)
            z = pose.p[2]

            #print desired_f[0, cnt]
            # account for the robot base and length of tool
            z = z + 0.026-0.18+0.00353
            v_z = Kp*(F-desired_f[0, cnt])
            v_l = np.array([v_x, v_y, v_z])

            force_out[:, cnt] = np.concatenate((T_b, F_b), axis=0)

            eef_pos[:, cnt] = pose.p

            quat = rox.R2q(R)
            eef_orien[:, cnt] = quat
            tim[0, cnt] = time.time()

            cnt += 1
     
        print F

        #### change here ####
        #pos = pos + v_l[2]*delta_t
        #acce = (v_l[2]-v_l_pre)/delta_t
						
        # formalize entire twist
        spatial_velocity_command = np.array([0, 0, 0, v_l[0], v_l[1], v_l[2]])

        # emergency stop if force too large
        if abs(F) > 2000:
            spatial_velocity_command = np.array([0, 0, 0, 0, 0, 0])
            print "force too large, stop..."
        
        # solve for joint velocity
        # Jacobian inverse
        #J = rox.robotjacobian(abb_robot, state.joint_angles)		
        #joints_vel = np.linalg.pinv(J).dot(spatial_velocity_command)
		
        # QP
        joints_vel = quadprog.compute_joint_vel_cmd_qp(state.joint_angles, spatial_velocity_command)
		
        # commanded joint position setpoint to EGM
        q_c = state.joint_angles + joints_vel*delta_t
		
        egm.send_to_robot(q_c)
   
        ####### change here ########
        #v_l_pre = v_l[2]
        
        #if t_new - t_pre < delta_t:
        #    time.sleep(delta_t - t_new + t_pre)	
		
        if cnt == n:
            csv_dat=np.hstack((desired_f.T, vdy.T, force_out.T, eef_pos.T, eef_orien.T, tim.T))
            np.savetxt('trap_force_trap_motion_020520.csv', csv_dat, fmt='%6.5f', delimiter=',')#, header='desired joint, optimal input')
            print "done"
				
		
if __name__ == '__main__':
    main()
