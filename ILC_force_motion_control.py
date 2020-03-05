import numpy as np
import rpi_abb_irc5
import time
import general_robotics_toolbox as rox
import collections
import scipy.io
import matplotlib.pyplot as plt
import QuadProg as qp
import sys
import rpi_ati_net_ft
from scipy.optimize import minimize_scalar
from scipy import interpolate
#import winsound

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
	
## does not track motion trajectory until reach the Fd
def first_half(input, num_iter):

    # stop the active RAPID program
    #rapid.stop()
    # reset the pointer to main
    #rapid.resetpp()
    print 'first half'
    print 'reset PP to main'
    time.sleep(5) 
    # start the RAPID program
    #rapid.start()
    
	# determine if robot has reached the initial configuration
    tag = True
    while tag:
        res, state = egm.receive_from_robot(.1)
        if res: 
            #print np.fabs(sum(np.rad2deg(state.joint_angles)) - sum(init))
            if np.fabs(sum(np.rad2deg(state.joint_angles)) - sum(init)) < 1e-4:
                tag = False
	
    time.sleep(1)	
    
    # out is composed of 5 velocity and 1 force in z
    out = np.zeros((6, n))
    force_out = np.zeros((6, n))		
    # pos of eef
    eef_pos = np.zeros((3, n))	
    # orientation of eef (quaternion)
    eef_orien = np.zeros((4, n))
    # timestamp
    tim = np.zeros((1, n))

    ############### change ################ or there will be errors that q_pre referred before assigned
    #q_pre = np.deg2rad(state.joint_angles)
    q_hat = np.zeros((6, 1))
    qhat_dot = np.zeros((6, 1))
    # for observer k should be symmetric and positive definite
    kl = 0.1

    ### drain the force sensor buffer ###
    count = 0
    while count < 1000:
        flag, ft = netft.read_ft_streaming(.1)
        #print ft[5]
        count = count+1     
    ### drain the EGM buffer ###
    for i in range(1000):
        res, state = egm.receive_from_robot(.1)
    # substract the initial force for bias
    flag, ft = netft.read_ft_streaming(.1)
    F0 = ft[5]
    print F0  
    time.sleep(3) 

    cnt = 0
    step_done = False
    while cnt < n:#pose.p[0] < 2 and cnt < n:
        #t_pre = time.time()
        # receive EGM feedback
        res, state = egm.receive_from_robot(.1)
		
        if not res:
            continue
			
        q_new = state.joint_angles

        if not step_done:
            print '--------start step-over motion--------'
            # do step-over of 0.25 mm in +x in world frame
            # current eef pose
            pose_cur = rox.fwdkin(abb_robot, q_new)
            pose_cur.p[0] = pose_cur.p[0] + num_iter*2*step_over
            # solve for inverse kinematics and pick the one that is closest to current configuration
            sol = rox.robot6_sphericalwrist_invkin(abb_robot, pose_cur, q_new)
            try:
                tar = sol[0] # raise exception if no solution
            except:
                tar = q_new
            # move to the new position after step-over
            egm.send_to_robot(tar)
            step_done = True

            q_new = tar
            ### drain the EGM buffer, or it will use the obsolete EGM feedback###
            for i in range(1000):
                res, state = egm.receive_from_robot(.1)
          
            print '--------step-over motion done--------'
            time.sleep(2)


	    # forward kinematics to calculate current position of eef
        pose = rox.fwdkin(abb_robot, q_new)
        R = pose.R
        flag, ft = netft.read_ft_streaming(.1)
        # map torque/force from sensor frame to base frame
        T_b = np.matmul(R, ft[0:3])
        F_b = np.matmul(R, ft[3:None])
        F = F_b[2]# - F0# first three torques and then three forces

        # start motion in y direction when robot barely touches coupon
        Fd0 = 50		
        if F < Fd0-0.1 and cnt == 0:
            z = pose.p[2]
            # account for the robot base and length of tool
            z = z + 0.026-0.18+0.00353
            # will shake if gain too large, here use 0.0002
            v_z = Kp*10*(F-Fd0)#-Ki*(z-z_ref)-Kd*acce#-Ki*pos
            # formalize entire twist
            spatial_velocity_command = np.array([0, 0, 0, 0, 0, v_z])
        else:
            # formalize entire twist
            # nominal input composed of F and v
            spatial_velocity_command = input[:, cnt] #np.array([0, 0, 0, vdx, 0, Fd])

            z = pose.p[2]
            # account for the robot base and length of tool
            z = z + 0.026-0.18+0.00353
            v_z = Kp*(F-spatial_velocity_command[5])#-Ki*(z-z_ref)-Kd*acce#-Ki*pos
            # nominal input only contains v
            spatial_velocity_command[5] = v_z
            # approximation of joint velocity
            #q_new = np.deg2rad(state.joint_angles)
			
            ######### change here, use observer instead of approximation to calculate q_dot ########
            if cnt == 0:
                q_hat = q_new
            qhat_dot = joints_vel + kl*(q_new-q_hat)
			
            #q_dot = (q_new - q_pre)/delta_t
            
            J = rox.robotjacobian(abb_robot, q_new)
            # estimate velocity
            v_est = J.dot(qhat_dot)
            # formalize the nominal output composed of F and v
            out[:, cnt] = np.append(v_est[0:5], F)
            			
            force_out[:, cnt] = np.concatenate((T_b, F_b), axis=0)

            eef_pos[:, cnt] = pose.p
            #eef_pos[2, cnt] = z
		
            R = pose.R
            quat = rox.R2q(R)
            eef_orien[:, cnt] = quat
            tim[0, cnt] = time.time()
            cnt = cnt+1

        print F
        # solve for joint velocity
        # Jacobian inverse
        #J = rox.robotjacobian(abb_robot, q_new)		
        #joints_vel = np.linalg.pinv(J).dot(spatial_velocity_command)

        # emergency stop if force too large
        if abs(F) > 2000:
            spatial_velocity_command = np.array([0, 0, 0, 0, 0, 0])
            print "force too large, stop..."

        # QP
        joints_vel = quadprog.compute_joint_vel_cmd_qp(q_new, spatial_velocity_command)
        
        # commanded joint position setpoint to EGM
        q_c = q_new + joints_vel*delta_t#*4
        egm.send_to_robot(q_c)
        # joint angle at previous time step
        #q_pre = q_new
        
        ############ change here ##############
        # make the time interval 0.004 s
        #t_new = time.time()
        #if t_new - t_pre < delta_t:
        #    time.sleep(delta_t - t_new + t_pre)	
           
		
        ######### change here ########
        q_hat = q_hat + qhat_dot*delta_t
			
    # interpolate to filter the output
    t_inter = np.arange(0, n*delta_t, delta_t)
    # interpolate each row of output
    for i in range(6):
        y = out[i, :]
        tck = interpolate.splrep(t_inter, y, s=0.01) # s = 0 means no interpolation
        ynew = interpolate.splev(t_inter, tck, der=0)
        out[i, :] = ynew
    
    error = out - desired
    # flip the error
    err_flip = np.fliplr(error)
    print np.linalg.norm(error, 'fro')
	
    return out, err_flip, np.linalg.norm(error, 'fro'), force_out, eef_pos, eef_orien, tim
    
def second_half(x, out_pre, num_iter):

    #rapid.stop()
    #rapid.resetpp()
    #time.sleep(2) 
    #print 'restart'
    #rapid.start()
    print 'second half'
    print 'reset PP to main'
    time.sleep(5) 
    # determine if robot has reached the initial configuration init
    tag = True
    while tag:
        res, state = egm.receive_from_robot(.1)
        if res: 
            #print np.fabs(sum(np.rad2deg(state.joint_angles)) - sum(init))
            if np.fabs(sum(np.rad2deg(state.joint_angles)) - sum(init)) < 1e-4:
                tag = False

    time.sleep(1)	
    
    out = np.zeros((6, n))

    ############### change ################ or there will be errors that q_pre referred before assigned
    #q_pre = np.deg2rad(state.joint_angles)
    q_hat = np.zeros((6, 1))
    qhat_dot = np.zeros((6, 1))
    # for observer k should be symmetric and positive definite
    kl = 0.1
	
    ### drain the force buffer ###
    count = 0
    while count < 1000:
        flag, ft = netft.read_ft_streaming(.1)
        #print ft[5]
        count = count+1        
    ### drain the EGM buffer ###
    for i in range(1000):
        res, state = egm.receive_from_robot(.1)

    # substract the initial force for bias
    flag, ft = netft.read_ft_streaming(.1)
    F0 = ft[5]
    print F0
    time.sleep(3)

    cnt = 0
    step_done = False
    while cnt < n:#pose.p[0] < 2 and cnt < n:
        #t_pre = time.time()
        # receive EGM feedback
        res, state = egm.receive_from_robot(.1)
        
        if not res:
            continue
       
        q_new = state.joint_angles
		
        if not step_done:
            print '--------start step-over motion--------'
            # do step-over of 0.25 mm in +x in world frame
            # current eef pose
            pose_cur = rox.fwdkin(abb_robot, q_new)
            pose_cur.p[0] = pose_cur.p[0] + (num_iter*2+1)*step_over
            # solve for inverse kinematics and pick the one that is closest to current configuration
            sol = rox.robot6_sphericalwrist_invkin(abb_robot, pose_cur, q_new)
            try:
                tar = sol[0] # raise exception if no solution
            except:
                tar = q_new
            # move to the new position after step-over
            egm.send_to_robot(tar)
            step_done = True

            q_new = tar
            ### drain the EGM buffer, or it will use the obsolete EGM feedback###
            for i in range(1000):
                res, state = egm.receive_from_robot(.1)

            print '--------step-over motion done--------'
            time.sleep(2)

	    # forward kinematics to calculate current position of eef
        pose = rox.fwdkin(abb_robot, q_new)
        R = pose.R
        flag, ft = netft.read_ft_streaming(.1)
        # map torque/force from sensor frame to base frame
        T_b = np.matmul(R, ft[0:3])
        F_b = np.matmul(R, ft[3:None])
        F = F_b[2]# - F0# first three torques and then three forces
		
        Fd0 = 50
        if F < Fd0-0.1 and cnt == 0:
            z = pose.p[2]
            # account for the robot base and length of tool
            z = z + 0.026-0.18+0.00353
            # will shake if gain too large, here use 0.0002
            v_z = Kp*10*(F-Fd0)#-Ki*(z-z_ref)-Kd*acce#-Ki*pos
            # formalize entire twist
            spatial_velocity_command = np.array([0, 0, 0, 0, 0, v_z])
        else:
            # formalize entire twist
            # nominal input composed of F and v
            spatial_velocity_command = x[:, cnt] #np.array([0, 0, 0, vdx, 0, Fd])
            z = pose.p[2]
            # account for the robot base and length of tool
            z = z + 0.026-0.18+0.00353
            v_z = Kp*(F-spatial_velocity_command[5])#-Ki*(z-z_ref)-Kd*acce#-Ki*pos
            # nominal input only contains v
            spatial_velocity_command[5] = v_z
            # approximation of joint velocity
            #q_new = np.deg2rad(state.joint_angles)
			
            ######### change here, use observer instead of approximation to calculate q_dot ########
            if cnt == 0:
                q_hat = q_new
            qhat_dot = joints_vel + kl*(q_new-q_hat)
							
            #q_dot = (q_new - q_pre)/delta_t
			
            J = rox.robotjacobian(abb_robot, q_new)
            # estimate velocity
            v_est = J.dot(qhat_dot)
            # formalize the nominal output composed of F and v
            out[:, cnt] = np.append(v_est[0:5], F) 
            cnt = cnt+1
		
        print F
        # solve for joint velocity
        # Jacobian inverse
        #J = rox.robotjacobian(abb_robot, q_new)		
        #joints_vel = np.linalg.pinv(J).dot(spatial_velocity_command)
		
        # emergency stop if force too large
        if abs(F) > 2000:
            spatial_velocity_command = np.array([0, 0, 0, 0, 0, 0])
            print "force too large, stop..."

        # QP
        joints_vel = quadprog.compute_joint_vel_cmd_qp(q_new, spatial_velocity_command)
			
        # commanded joint position setpoint to EGM
        q_c = q_new + joints_vel*delta_t#*4
		
        egm.send_to_robot(q_c)
        # joint angle at previous time step
        #q_pre = q_new
        #t_pre = t_new
        # make the time interval 0.004 s
        #t_new = time.time()
        #if t_new - t_pre < delta_t:
        #    time.sleep(delta_t - t_new + t_pre)	

        ######### change here ########
        q_hat = q_hat + qhat_dot*delta_t
		
		
    # interpolate to filter the output
    t_inter = np.arange(0, n*delta_t, delta_t)
    # interpolate each row of output
    for i in range(6):
        y = out[i, :]
        tck = interpolate.splrep(t_inter, y, s=0.01) # s = 0 means no interpolation
        ynew = interpolate.splev(t_inter, tck, der=0)
        out[i, :] = ynew
        
    err = out-out_pre
    err_flip2 = np.fliplr(err)
    
    return err_flip2


# initialize EGM interface instance
egm = rpi_abb_irc5.EGM()

# initialize a robot instance
abb_robot = abb_irb6640_180_255_robot()

# desired force
#Fd = -500

# desired velocity in y
#vdy = -0.5

# feedback gain
Kp = 0.0002
Kd = 0.0008
Ki = 0.0004

# time step
delta_t = 0.004

# initial configuration in degree
init = [-91.08,2.54,38.18,0.0,49.27,-1.07]

step_over = 0.00025 # in meter

# quadprog to solve for joint velocity
quadprog = qp.QuadProg(abb_robot)

try:
    if (len(sys.argv) < 2):
        raise Exception('IP address of ATI Net F/T sensor required')
    host=sys.argv[1]
    netft=rpi_ati_net_ft.NET_FT(host)
    netft.set_tare_from_ft()
   
    netft.start_streaming()
    
except KeyboardInterrupt:
    pass

pos = 0
acce = 0
v_l_pre = 0

############ change here, how long the process lasts ############
n = 2000

####### trapezoidal desired force in z #######
tt = np.linspace(0, 4*np.pi, n)
Fdz = np.zeros((n, ))
vdy = np.zeros((n, ))
# form trap force and trap motion
for i in range(n):
    if tt[i] >= 0 and tt[i] < np.pi:
        Fdz[i] = 50+302*tt[i]
        vdy[i] = -0.2*tt[i]
    elif tt[i] >= np.pi and tt[i] < 3*np.pi:
        Fdz[i] = 50+302*np.pi
        vdy[i] = -0.2*np.pi
    else:
        Fdz[i] = 50+302*np.pi-302*(tt[i]-3*np.pi)
        vdy[i] = -0.2*np.pi + 0.2*(tt[i]-3*np.pi)

#Fdz = np.zeros((n, )) 
#plt.plot(tt, Fdz)
#plt.show()

#vdy = 0.5*2*np.sin(5*tt)

# received output (composed of force and velocity)
#out = np.zeros((6, n))

# desired output (composed of force and velocity)
desired = np.zeros((6, n))

# form the desired output by ud
for i in range(n):
    ud = np.array([0, 0, 0, 0, vdy[i], Fdz[i]])
    desired[:, i] = ud

# referece height of coupon that achieves desired force
z_ref = 0.89226
x_ref = 0
y_ref = -1.35626
	
# RobotStudio network service
#if (len(sys.argv) >= 2):
#    rapid=rpi_abb_irc5.RAPID(sys.argv[2])
#else:
#    rapid=rpi_abb_irc5.RAPID()

# input to the nominal dynamical system that is subject to update by ILC
# at the beginning, it is just desired output
desired_cp = desired.copy()
x_in = desired_cp

fro_err_old = 0

####### change here #######
iter = 20
for i in range(iter):
    # first pass into the dynamical system
    x_in_cp = x_in.copy()
    out, err_flip1, fro_err, force_out, eef_pos, eef_orien, tim = first_half(x_in_cp, i)
    
    # save all data after each iteration
    csv_dat=np.hstack((desired.T, x_in.T, out.T, force_out.T, eef_pos.T, eef_orien.T, tim.T))
    np.savetxt('ILC_trap_force_trap_motion_control_with_' + str(i) + '_iteration.csv', csv_dat, fmt='%6.5f', delimiter=',')
     
    #print "done"
    # check if the stopping condition satisfied
    #if np.fabs(fro_err-fro_err_old) < 0.5 or i == iter-1:
    #    csv_dat=np.hstack((desired.T, x_in.T, out.T))
		
        ########################
        ######change here#######
    #    np.savetxt('force_motion_control_final_sin_motion_101419.csv', csv_dat, fmt='%6.5f', delimiter=',')#, header='desired joint, optimal input')
    #    print "done"
        #frequency = 2500  # Set Frequency To 2500 Hertz
        #duration = 1000  # Set Duration To 1000 ms == 1 second
        #winsound.Beep(frequency, duration)
    #    break;
	
    time.sleep(2)   
    
    x = x_in+err_flip1	
    x_cp = x.copy()
    # second pass into the dynamical system	
    errflip2 = second_half(x_cp, out, i)
    #plt.plot(tt, errflip2[5, :])
    #plt.show()
    time.sleep(2)  
    
	################### sometimes 1d search may still be slow
    #print '----start searching optimal learning rate.----'
    #res = minimize_scalar(object_function, bounds=(0.0, 1.0), method='bounded', options={'maxiter': 5})
    #print '----the optimal learning rate is:----' 
    #print res.x
    
    # update input
    #x_in = x_in - res.x*errflip2
	
    ############### use fixed learning rate
    res = 0.25
    x_in = x_in - res*errflip2
	
    fro_err_old = fro_err
