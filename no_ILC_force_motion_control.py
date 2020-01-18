import numpy as np
import rpi_abb_irc5
import time
import general_robotics_toolbox as rox
import QuadProg as qp
import sys
import rpi_ati_net_ft

def abb_irb6640_180_255_robot():
    """Return a Robot instance for the ABB IRB6640 180-255 robot"""
    
    x = np.array([1,0,0])
    y = np.array([0,1,0])
    z = np.array([0,0,1])
    a = np.array([0,0,0])
    
    H = np.array([z,y,y,x,y,x]).T
    P = np.array([0.78*z, 0.32*x, 1.075*z, 0.2*z, 1.142*x, 0.2*x, a]).T
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
    Fd = -1000
		
    # desired velocity in y
    vdy = 0
	
    # feedback gain
    Kp = 0.0004
    Kd = 0.0008*2
    Ki = 0.0004
	
    # time step
    delta_t = 0.004
	
    # initial configuration in degree
    init = [-90,1.92,38.8,0,48.28,0]

    # quadprog to solve for joint velocity
    quadprog = qp.QuadProg(abb_robot)
	
    # force sensor initialization
    try:
        if (len(sys.argv) < 2):
            raise Exception('IP address of ATI Net F/T sensor required')
        host=sys.argv[1]
        netft=rpi_ati_net_ft.NET_FT(host)
        netft.set_tare_from_ft()
        
        netft.start_streaming()
       
    except KeyboardInterrupt:
        pass
	
    ######## change here ########
    pos = 0
    acce = 0
    v_l_pre = 0
    #########
	
    # output force in z
    n = 6000
    force_out = np.zeros((6, n))	
	
    # desired force in z
    desired_f = np.zeros((1, n))
	
    # pos of eef
    eef_pos = np.zeros((3, n))	
	
    # referece height of coupon that achieves desired force
    z_ref = 0.89226
    x_ref = 0
    y_ref = -1.35626
	
    # determine if robot has reached the initial configuration init
    tag = True
    while tag:
        res, state = egm.receive_from_robot(.1)
        if res: 
            if np.fabs(sum(np.rad2deg(state.joint_angles)) - sum(init)) < 1e-2:
                tag = False
	
    time.sleep(0.5)	
    print '--------start force control--------'
	
    cnt = 0
    count = 0
    while count < n: #pose.p[0] < 2:
        t_pre = time.time()
        # receive EGM feedback
        res, state = egm.receive_from_robot(.1)
        #print state.joint_angles
        if not res:
            continue
       
	    # forward kinematics to calculate current position of eef
        pose = rox.fwdkin(abb_robot, state.joint_angles)
        # read force/torque
        flag, ft = netft.read_ft_streaming(.1)
        F = ft[5] # first three torques and then three forces

        if F > Fd+0.1 and cnt == 0:
            Fd = Fd
            z = pose.p[2]
            # account for the robot base and length of tool
            z = z + 0.026-0.18
            # will shake if gain too large, here use 0.0002
            v_z = -0.0002*1*(F-Fd)-Ki*(z-z_ref)-Kd*acce#-Ki*pos
            v_l = np.array([0, 0, v_z])
        else: 
            Fd = Fd
            # ramping up force

            cnt += 1
            if cnt == 1:
                tt = np.linspace(0, 4*np.pi, n)
                vdy = 0.5*2*np.sin(5*tt)
                #print "start motion"
            
            z = pose.p[2]
            # account for the robot base and length of tool
            z = z + 0.026-0.18
            v_z = -Kp*(F-Fd)-Ki*(z-z_ref)-Kd*acce
            #v_z = -Kp*1.0*(F-Fd)-Kd*acce#-Ki*pos
            v_l = np.array([0, vdy[cnt-1], v_z])
          
        print F
        force_out[:, count] = ft
        eef_pos[:, count] = pose.p
        eef_pos[2, count] = pose.p[2] + 0.026-0.18
        desired_f[0, count] = Fd
		
        #### change here ####
        #pos = pos + v_l[2]*delta_t
        acce = (v_l[2]-v_l_pre)/delta_t
						
        # formalize entire twist
        spatial_velocity_command = np.array([0, 0, 0, v_l[0], v_l[1], v_l[2]])

        # emergency stop if force too large
        if abs(F) > 2500:
            spatial_velocity_command = np.array([0, 0, 0, 0, 0, 0])
            print "force too large, stop..."
        # solve for joint velocity
        # Jacobian inverse
        #J = rox.robotjacobian(abb_robot, np.deg2rad(state.joint_angles))		
        #joints_vel = np.linalg.pinv(J).dot(spatial_velocity_command)
		
        # QP
        joints_vel = quadprog.compute_joint_vel_cmd_qp(state.joint_angles, spatial_velocity_command)
		
        # commanded joint position setpoint to EGM
        q_c = state.joint_angles + joints_vel*delta_t*4
		
        egm.send_to_robot(q_c)
   
        ####### change here ########
        v_l_pre = v_l[2]
        
        t_new = time.time()

        #t_all = t_all+ t_new-t_pre
        if t_new - t_pre < delta_t:
            time.sleep(delta_t - t_new + t_pre)	
        count = count+1
		
        if count == n:
            csv_dat=np.hstack((desired_f.T, force_out.T, eef_pos.T))
            np.savetxt('constant_force_slow_motion_mass_16_120419_same_KpKi.csv', csv_dat, fmt='%6.5f', delimiter=',')#, header='desired joint, optimal input')
            print "done"
				
		
if __name__ == '__main__':
    main()
