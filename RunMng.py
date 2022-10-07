from scipy.linalg import solve
import time
from scipy import integrate
import scipy as sp
import scipy
import machine_constant_config as mc
import sys,os
sys.path.append(os.path.abspath(r'E:\PycharmProject\pythonProject'))
from mypylib import *
import control as ctrl

gl_servo_freq=0
gl_T = 0
gl_plant=0

"""
Plant Init
"""
def init_plant():
    # 建立位置环受控对象
    global gl_plant
    m = mc.plant.mass
    c = mc.plant.c
    k = mc.plant.k
    A = np.array([[0, 1], [-k / m, -c / m]])
    B = np.array([[0], [1 / m]])
    C = np.array([1, 0])
    D = 0.0
    gl_plant = ctrl.StateSpace(A, B, C, D)

def init_plant_with_resonance():
    global gl_plant
    m = mc.plant.mass

    zd=0.1
    zf=200
    pd=0.1
    pf=300
    zw = zf * 2 * np.pi
    pw = pf * 2 * np.pi
    k = pw**2 / zw**2
    num_s = k * np.array([1, 2 * zd * zw, zw**2])
    den_s = np.array([1, 2 * pd * pw, pw**2])
    den_s=np.append(den_s*m,[0,0])
    A, B, C, D = signal.tf2ss(num_s, den_s)
    # A, B, C, D = ctrl.tf2ss(num_s, den_s)
    gl_plant = ctrl.StateSpace(A, B, C, D)

    tf=ctrl.tf(num_s,den_s)
    w=np.arange(30,1000,1/10000)*2*np.pi
    mag, phase, omega =ctrl.bode_plot(tf,w,dB=True,Hz=True,deg=True,Plot=True)



"""
PID Init
"""
def pid_design(method,mass,dampling=0.707,ratio=0.5,alpha=0.5,bandwidth_Hz=100):
    if method == "pole_placement":
        omega = 2*np.pi*bandwidth_Hz
        k = 1 / mass
        kp = (1 + 2 * dampling * ratio) * omega ** 2 / k
        ki = ratio * omega ** 3 / k / gl_servo_freq#一阶后向系数
        kd = (2 * dampling * omega + ratio * omega - alpha) / k * gl_servo_freq#一阶后向系数
        # ki = ratio * omega ** 3 / k  * gl_servo_freq*2#双边线性系数
        # kd = (2 * dampling * omega + ratio * omega - alpha) / k / gl_servo_freq/2#双边线性系数

    elif method == "hans_buttler":
        alpha = 3
        wn_bw = 2 * np.pi * bandwidth_Hz
        omega_d = wn_bw / alpha
        omega_i = wn_bw / alpha ** 2
        kp_series = mass * wn_bw ** 2 / alpha
        kp = kp_series * (omega_i / omega_d + 1)
        ki = kp_series * omega_i / gl_servo_freq#一阶后向系数
        kd = kp_series / omega_d * gl_servo_freq#一阶后向系数

        # ki = kp_series * omega_i * gl_servo_freq*2#双边线性系数
        # kd = kp_series / omega_d / gl_servo_freq/2#双边线性系数
    return kp, ki, kd

"""
SPG Init
"""
def gen_spg(method,time,distance,settling_time):
    if method == "5_order":
        t_list = np.linspace(0, time, int(time * gl_servo_freq))
        #p0=0,v0=0,a0=0, v1=0,a1=0
        p1,v1,a1,t1=distance,0,0,time
        A=np.array([[t1**5,t1**4,t1**3],
                   [5*t1**4,4*t1**3,3*t1**2],
                   [20*t1**3,12*t1**2,6*t1]])
        b=np.array([p1,v1,a1])
        X=solve(A,b)
        set_pos=0.0+X[0]*t_list**5+X[1]*t_list**4+X[2]*t_list**3
        set_vel=0.0+5*X[0]*t_list**4+4*X[1]*t_list**3+3*X[2]*t_list**2
        set_acc=0.0+20*X[0]*t_list**3+12*X[1]*t_list**2+6*X[2]*t_list
        set_jerk=0.0+60*X[0]*t_list**2+24*X[1]*t_list+6*X[2]
    elif method=="sin":
        pass
    settling_time_count=int(settling_time*gl_servo_freq)
    set_pos=np.pad(set_pos, (0, settling_time_count), 'edge')
    set_vel=np.pad(set_vel, (0, settling_time_count), 'constant')
    set_acc=np.pad(set_acc, (0, settling_time_count), 'constant')
    set_jerk=np.pad(set_jerk, (0, settling_time_count), 'constant')
    return set_pos,set_vel,set_acc,set_jerk

def run_closed_loop(enable_save_datalog=False):
    plt.close("all")
    # 初始化控制器
    # Kp, Ki, Kd = pid_design(method="pole_placement", mass=mc.plant.mass, dampling=0.707, ratio=0.5, alpha=mc.plant.c/mc.plant.mass, bandwidth_Hz=100)
    Kp, Ki, Kd = pid_design(method="hans_buttler", mass=mc.plant.mass, dampling=0.707, ratio=0.5,
                            alpha=mc.plant.c / mc.plant.mass, bandwidth_Hz=50)
    # 初始化SPG
    profile, set_vel, set_acc, set_jerk = gen_spg(method="5_order", time=1, distance=0.2,settling_time=0.1)
    total_samples = len(profile)
    pos_fbk = 0  # 当前起始位置
    enc_fbk = 0  # 当前起始位置
    vel_fbk = 0  # 当前起始位置
    last_ctrl_output = 0
    interateout = 0  #
    differential = 0  #
    err = np.zeros(2)
    Datalog = np.zeros((total_samples, 10))
    # begin runmng
    for i in range(total_samples):
        err[1] = profile[i] - enc_fbk
        Datalog[i, 1] = err[1]
        interateout += Ki * err[1]#一阶后向差分
        # interateout += Ki*(err[0]+err[1])#双边线性差分
        # 积分限幅
        interateout = 1950 if interateout > 1950 else interateout
        interateout = -1950 if interateout < -1950 else interateout

        differential=Kd * (err[1] - err[0])#一阶后向差分
        # differential=Kd*(err[1] - err[0])-differential#双边线性差分

        ctrlOutput = Kp * err[1] + interateout + differential
        # 控制器输出限幅
        ctrlOutput = 2000 if ctrlOutput > 2000 else ctrlOutput
        ctrlOutput = -2000 if ctrlOutput < -2000 else ctrlOutput
        err[0] = err[1]
        # 仿真实物运动
        # tout, y, x = ctrl.forced_response(plant, T=[0, gl_T], U=[ctrlOutput, ctrlOutput], X0=[enc_fbk, vel_fbk],return_x=True)  # 使用的是zoh法,control库版本 0.9.2
        tout, y, x = ctrl.forced_response(gl_plant, T=[0, gl_T], U=[ctrlOutput, ctrlOutput],X0=[enc_fbk, vel_fbk])  # control库0.8.3
        vel_fbk = x[-1][-1]
        enc_fbk = y[-1]
        last_ctrl_output = ctrlOutput

        Datalog[i, 0] = i
        Datalog[i, 3] = enc_fbk
        Datalog[i, 2] = profile[i]
        Datalog[i, 4] = ctrlOutput
        Datalog[i, 6] = set_vel[i]
        Datalog[i, 7] = set_acc[i]
        Datalog[i, 8] = vel_fbk
        Datalog[i, 9] = set_jerk[i]

    # t, series_pos_fbk,xout = ctrl.forced_response(gl_plant, T=Datalog[:, 0] * gl_T, U=Datalog[:, 4], X0=0)
    # Datalog[:, 5] = series_pos_fbk
    # 绘制结果
    plot_uy(Datalog[:, 2], Datalog[:, 3], "set_pos", "enc", 2, 1)
    plot_uy(Datalog[:, 4], Datalog[:, 1], "ctrlOutput", "pos_err", 2, 1)
    # 存储结果
    if enable_save_datalog:
        now_time = time.strftime('%Y_%m_%d_%H_%M_%S', time.localtime())
        path = os.getcwd()
        file_path = path + "\\Datalog\\" + now_time + '.csv'
        df1 = pd.DataFrame(data=Datalog,
                           columns=['sample', 'pos_err', 'set_pos', 'enc', 'ctrlOutput', 'series_pos_fbk', 'set_vel',
                                    'set_acc', 'vel_fbk', 'set_jerk'])
        df1.to_csv(file_path, index=False)


def set_machine_constant():
    global gl_servo_freq,gl_T
    gl_servo_freq = mc.system.gl_servo_freq
    gl_T = 1 / gl_servo_freq
    init_plant()
    # init_plant_with_resonance()
def chirp_signal_plant_open_loop_acc(enable_save_datalog=False):
    start_freq = 30
    stop_freq = 4000
    Force=40
    T = 2#时间s

    f_chirp_start=start_freq*0.8
    f_chirp_end=stop_freq*1.2
    duration=(f_chirp_end-f_chirp_start)/(stop_freq-start_freq)*T
    t_list =np.arange(0,duration,1/gl_servo_freq)#step =1/servo_freq_current

    F_disturbance = 0#wgn(F_disturbance, snr=10)  # 给构造电流采集信号加白噪声,snr信噪比

    Datalog = np.zeros((len(t_list), 4))
    chirp_cmd = Force*np.sin(2*np.pi*((f_chirp_end-f_chirp_start)/duration*t_list**2/2+f_chirp_start*t_list))
    # force=F_disturbance+chirp_cmd
    # foece_noise=wgn(chirp_cmd, snr=30)
    # t, pos_fbk, xout = ctrl.forced_response(gl_plant, T=t_list, U=foece_noise, X0=0)

    input=chirp_cmd#wgn(chirp_cmd, snr=30)
    # begin runmng
    t, yout, xout = ctrl.forced_response(gl_plant, T=t_list, U=input, X0=0)
    output=yout#wgn(yout, snr=30)  # 给构造传感器白噪声,snr信噪比
    plot_uy(input, output, "chirp_cmd", "yout")
    acc=np.diff(output,2)*gl_servo_freq*gl_servo_freq
    u=input
    y=np.pad(acc,(2,0),"constant",constant_values=0)
    plot_uy(input, y, "chirp_cmd", "acc")
    fw,acc_fw=sys_id_correlate_resolution(gl_servo_freq, u, y, f_chirp_start/0.8, f_chirp_end/1.2,resolution=0.5)
    # fw,acc_fw=sys_id_fft(gl_servo_freq, u, y, f_chirp_start/0.8, f_chirp_end/1.2)
    s_d = 1j * 2 * np.pi * fw
    decay=(1-np.e**(-s_d*gl_T))/s_d/gl_T
    acc_fw=acc_fw/decay/decay
    pos_fw=acc_fw/s_d/s_d
    # bode_plot(fw,acc_fw,"Acc")
    bode_plot(fw,pos_fw,"Pos")
    # plot_uy(y, foece_noise, "acc", "foece_noise", 2, 1)
    # plot_uy(pos_fbk, pos_fbk_noise, "pos_fbk", "pos_fbk_noise", 2, 1)
    # plot_uy(chirp_cmd, F_disturbance, "chirp_cmd", "F_disturbance", 2, 1)

    # 存储结果
    if enable_save_datalog:
        now_time = time.strftime('%Y_%m_%d_%H_%M_%S', time.localtime())
        path = os.getcwd()
        file_path = path + "\\Datalog\\" + now_time + '.csv'
        df1 = pd.DataFrame(data=Datalog,
                           columns=['sample', 'chirp_cmd', 'acc','pos_fbk'])
        df1.to_csv(file_path, index=False)

def chirp_signal_plant_open_loop(enable_save_datalog=False):
    start_freq = 30
    stop_freq = 1000
    Force=40
    T = 2#时间s
    f_chirp_start=start_freq*0.8
    f_chirp_end=stop_freq*1.2
    duration=(f_chirp_end-f_chirp_start)/(stop_freq-start_freq)*T
    t_list =np.arange(0,duration,1/gl_servo_freq)#step =1/servo_freq_current
    chirp_cmd = Force*np.sin(2*np.pi*((f_chirp_end-f_chirp_start)/duration*t_list**2/2+f_chirp_start*t_list))
    total_samples=len(t_list)
    input=chirp_cmd#wgn(chirp_cmd, snr=30)
    # begin runmng
    t, yout, xout = ctrl.forced_response(gl_plant, T=t_list, U=input, X0=0)
    output=yout#wgn(yout, snr=30)  # 给构造传感器白噪声,snr信噪比
    u=input
    y=output
    plot_uy(u, y, "input", "output")
    # fw,transfer_fw=sys_id_correlate_resolution(gl_servo_freq, u, y, f_chirp_start/0.8, f_chirp_end/1.2,resolution=0.5)
    fw,transfer_fw=sys_id_fft(gl_servo_freq, u, y, f_chirp_start/0.8, f_chirp_end/1.2)
    bode_plot(fw,transfer_fw,"transfer_fw")

def chirp_signal_plant_after_controller_close_loop(enable_save_datalog=False):
    plt.close("all")
    # 初始化扫频参数
    start_freq = 10
    stop_freq = 4000
    chirp_force=40
    T = 1#时间s
    f_chirp_start=start_freq*0.8
    f_chirp_end=stop_freq*1.1
    duration=(f_chirp_end-f_chirp_start)/(stop_freq-start_freq)*T
    t_list =np.arange(0,duration,1/gl_servo_freq)#step =1/servo_freq_current
    total_samples=len(t_list)
    chirp_cmd = chirp_force*np.sin(2*np.pi*((f_chirp_end-f_chirp_start)/duration*t_list**2/2+f_chirp_start*t_list))
    # 初始化控制器
    Kp, Ki, Kd = pid_design(method="hans_buttler", mass=mc.plant.mass, dampling=0.707, ratio=0.5,
                            alpha=mc.plant.c / mc.plant.mass, bandwidth_Hz=50)
    interateout = 0  #
    enc_fbk = 0  # 当前起始位置
    vel_fbk = 0  # 当前起始位置
    err = np.zeros(2)
    Datalog = np.zeros((total_samples, 10))
    profile=np.zeros((total_samples, 1))
    ffc_out=0
    F_disturbance=0#初始化外界干扰力
    # begin runmng
    for i in range(total_samples):
        err[1] = profile[i] - enc_fbk
        Datalog[i, 1] = err[1]
        interateout += Ki * err[1]
        # 积分限幅
        interateout = 1950 if interateout > 1950 else interateout
        interateout = -1950 if interateout < -1950 else interateout

        pid_put = Kp * err[1] + interateout + Kd * (err[1] - err[0])
        # 控制器输出限幅
        pid_put = 2000 if pid_put > 2000 else pid_put
        pid_put = -2000 if pid_put < -2000 else pid_put
        err[0] = err[1]
        ctrl_out=pid_put+ffc_out+chirp_cmd[i]


        # 仿真实物运动
        # F_disturbance=10**5*enc_fbk+20  #模拟SSZ向上磁弹力
        force=ctrl_out+F_disturbance
        tout, y, x = ctrl.forced_response(gl_plant, T=[0, gl_T], U=[force, force],
                                          X0=[enc_fbk, vel_fbk])  # control库0.8.3
        vel_fbk = x[-1][-1]
        enc_fbk = y[-1]

        Datalog[i, 0] = i
        Datalog[i, 3] = enc_fbk
        Datalog[i, 2] = profile[i]
        Datalog[i, 4] = force
        Datalog[i, 5] = ctrl_out
        Datalog[i, 6] = chirp_cmd[i]
        Datalog[i, 7] = pid_put
        Datalog[i, 8] = vel_fbk
        Datalog[i, 9] = F_disturbance
    # 绘制结果
    plot_uy(Datalog[:, 2], Datalog[:, 3], "set_pos", "enc", 2, 1)
    plot_uy(Datalog[:, 4], Datalog[:, 1], "force", "pos_err", 2, 1)
    u=Datalog[:, 5]+Datalog[:, 9]#𝑃𝐼𝐷+𝑆𝑖𝑔𝑖𝑛+𝐷𝑖𝑠𝑡𝑢𝑟𝑏
    y=-Datalog[:, 7]#-𝑃𝐼𝐷
    plot_uy(u, y, "ctrl_out+Disturb", "-PID")#所求为CP的传函
    plot_uy(Datalog[:, 5], Datalog[:, 9], "ctrl_out", "Disturb")#所求为CP的传函
    fw,transfer_fw=sys_id_correlate_resolution(gl_servo_freq, u, y, f_chirp_start/0.8, f_chirp_end/1.1,resolution=1.5)
    # fw,acc_fw=sys_id_fft(gl_servo_freq, u, y, f_chirp_start, f_chirp_end)
    s_d = 1j * 2 * np.pi * fw
    z=np.e**(s_d*gl_T)
    s_z=1-z**(-1)#一阶后向差分
    c_fw=Kp+Ki*1/s_z+Kd*s_z
    pos_fw=transfer_fw/c_fw
    close_loop_fw=transfer_fw/1+transfer_fw#transfer_fw所求为CP的传函
    acc_fw=pos_fw*s_d*s_d
    bode_plot(fw,c_fw,"c_fw")
    bode_plot(fw,transfer_fw,"cp")
    bode_plot(fw,close_loop_fw,"close_loop_fw")
    bode_plot(fw,acc_fw,"Acc")
    bode_plot(fw,pos_fw,"Pos")

    # 存储结果
    if enable_save_datalog:
        now_time = time.strftime('%Y_%m_%d_%H_%M_%S', time.localtime())
        path = os.getcwd()
        file_path = path + "\\Datalog\\" + now_time + '.csv'
        df1 = pd.DataFrame(data=Datalog,
                           columns=['sample', 'pos_err', 'set_pos', 'enc', 'force', 'ctrl_out', 'chirp_cmd',
                                    'pid_out', 'vel_fbk', 'set_jerk'])
        df1.to_csv(file_path, index=False)

def chirp_signal_plant_after_controller_close_loop_sigin_pidout(enable_save_datalog=False):
    plt.close("all")
    # 初始化扫频参数
    start_freq = 30
    stop_freq = 4000
    chirp_force=40
    T = 1#时间s
    f_chirp_start=start_freq*0.8
    f_chirp_end=stop_freq*1.2
    duration=(f_chirp_end-f_chirp_start)/(stop_freq-start_freq)*T
    t_list =np.arange(0,duration,1/gl_servo_freq)#step =1/servo_freq_current
    total_samples=len(t_list)
    chirp_cmd = chirp_force*np.sin(2*np.pi*((f_chirp_end-f_chirp_start)/duration*t_list**2/2+f_chirp_start*t_list))
    # 初始化控制器
    Kp, Ki, Kd = pid_design(method="hans_buttler", mass=mc.plant.mass, dampling=0.707, ratio=0.5,
                            alpha=mc.plant.c / mc.plant.mass, bandwidth_Hz=50)
    interateout = 0  #
    enc_fbk = 0  # 当前起始位置
    vel_fbk = 0  # 当前起始位置
    err = np.zeros(2)
    Datalog = np.zeros((total_samples, 10))
    profile=np.zeros((total_samples, 1))
    ffc_out=0
    F_disturbance=0#初始化外界干扰力
    # begin runmng
    for i in range(total_samples):
        err[1] = profile[i] - enc_fbk
        Datalog[i, 1] = err[1]
        interateout += Ki * err[1]
        # 积分限幅
        interateout = 1950 if interateout > 1950 else interateout
        interateout = -1950 if interateout < -1950 else interateout

        pid_put = Kp * err[1] + interateout + Kd * (err[1] - err[0])
        # 控制器输出限幅
        pid_put = 2000 if pid_put > 2000 else pid_put
        pid_put = -2000 if pid_put < -2000 else pid_put
        err[0] = err[1]
        ctrl_out=pid_put+ffc_out+chirp_cmd[i]
        # 仿真实物运动
        # F_disturbance=25  #模拟SSZ向上磁弹力
        force=ctrl_out+F_disturbance
        tout, y, x = ctrl.forced_response(gl_plant, T=[0, gl_T], U=[force, force],
                                          X0=[enc_fbk, vel_fbk])  # control库0.8.3
        vel_fbk = x[-1][-1]
        enc_fbk = y[-1]

        Datalog[i, 0] = i
        Datalog[i, 2] = profile[i]
        Datalog[i, 3] = enc_fbk
        Datalog[i, 4] = force
        Datalog[i, 5] = ctrl_out
        Datalog[i, 6] = chirp_cmd[i]
        Datalog[i, 7] = pid_put
        Datalog[i, 8] = vel_fbk
        Datalog[i, 9] = F_disturbance
    # 绘制结果
    plot_uy(Datalog[:, 2], Datalog[:, 3], "set_pos", "enc", 2, 1)
    plot_uy(Datalog[:, 4], Datalog[:, 1], "force", "pos_err", 2, 1)
    u=chirp_cmd
    y=Datalog[:, 7]
    # u = u-np.mean(u)
    # y = y-np.mean(y)
    plot_uy(u, y, "chirp_cmd", "pid_put")#
    fw,transfer_fw=sys_id_correlate_resolution(gl_servo_freq, u, y, f_chirp_start/0.8, f_chirp_end/1.2,resolution=1.5)
    # fw,transfer_fw=sys_id_fft(gl_servo_freq, u, y, f_chirp_start/0.8, f_chirp_end/1.2)
    s_d = 1j * 2 * np.pi * fw
    z=np.e**(s_d*gl_T)
    s_z=1-z**(-1)#一阶后向差分,此处不除以T是因为Ki和Kd参数已经除过了
    c_fw=Kp+Ki*1/s_z+Kd*s_z
    close_loop_fw=-transfer_fw
    cp=close_loop_fw/(np.ones_like(fw)-close_loop_fw)
    pos_fw=cp/c_fw
    acc_fw=pos_fw*s_d*s_d
    bode_plot(fw,c_fw,"c_fw")
    bode_plot(fw,cp,"cp")
    bode_plot(fw,close_loop_fw,"close_loop_fw")
    bode_plot(fw,acc_fw,"Acc")
    bode_plot(fw,pos_fw,"Pos")

    # 存储结果
    if enable_save_datalog:
        now_time = time.strftime('%Y_%m_%d_%H_%M_%S', time.localtime())
        path = os.getcwd()
        file_path = path + "\\Datalog\\" + now_time + '.csv'
        df1 = pd.DataFrame(data=Datalog,
                           columns=['sample', 'pos_err', 'set_pos', 'enc', 'force', 'ctrl_out', 'chirp_cmd',
                                    'pid_out', 'vel_fbk', 'set_jerk'])
        df1.to_csv(file_path, index=False)
def chirp_signal_plant_after_controller_close_loop_sigin_ctrlout(enable_save_datalog=False):
    plt.close("all")
    # 初始化扫频参数
    start_freq = 10
    stop_freq = 2000
    chirp_force=40
    T = 2#时间s
    f_chirp_start=start_freq*0.8
    f_chirp_end=stop_freq*1.2
    duration=(f_chirp_end-f_chirp_start)/(stop_freq-start_freq)*T
    t_list =np.arange(0,duration,1/gl_servo_freq)#step =1/servo_freq_current
    total_samples=len(t_list)
    chirp_cmd = chirp_force*np.sin(2*np.pi*((f_chirp_end-f_chirp_start)/duration*t_list**2/2+f_chirp_start*t_list))
    # 初始化控制器
    Kp, Ki, Kd = pid_design(method="hans_buttler", mass=mc.plant.mass, dampling=0.707, ratio=0.5,
                            alpha=mc.plant.c / mc.plant.mass, bandwidth_Hz=50)
    interateout = -19.54  #
    enc_fbk = 0  # 当前起始位置
    vel_fbk = 0  # 当前起始位置
    err = np.zeros(2)
    Datalog = np.zeros((total_samples, 11))
    profile=np.zeros((total_samples, 1))
    ffc_out=0
    F_disturbance=-1*interateout#初始化外界干扰力
    # begin runmng
    for i in range(total_samples):
        err[1] = profile[i] - enc_fbk
        Datalog[i, 1] = err[1]
        interateout += Ki * err[1]
        # 积分限幅
        interateout = 1950 if interateout > 1950 else interateout
        interateout = -1950 if interateout < -1950 else interateout

        pid_put = Kp * err[1] + interateout + Kd * (err[1] - err[0])
        # 控制器输出限幅
        pid_put = 2000 if pid_put > 2000 else pid_put
        pid_put = -2000 if pid_put < -2000 else pid_put
        err[0] = err[1]
        ctrl_out=pid_put+ffc_out+chirp_cmd[i]
        # 仿真实物运动
        F_disturbance=enc_fbk+19.582  #模拟SSZ向上磁弹力
        force=ctrl_out+F_disturbance
        tout, y, x = ctrl.forced_response(gl_plant, T=[0, gl_T], U=[force, force],
                                          X0=[enc_fbk, vel_fbk])  # control库0.8.3
        vel_fbk = x[-1][-1]
        enc_fbk = y[-1]

        Datalog[i, 0] = i
        Datalog[i, 2] = profile[i]
        Datalog[i, 3] = enc_fbk
        Datalog[i, 4] = force
        Datalog[i, 5] = ctrl_out
        Datalog[i, 6] = chirp_cmd[i]
        Datalog[i, 7] = pid_put
        Datalog[i, 8] = vel_fbk
        Datalog[i, 9] = F_disturbance
        Datalog[i, 10] = interateout
    # 绘制结果
    plot_uy(Datalog[:, 2], Datalog[:, 3], "set_pos", "enc", 2, 1)
    plot_uy(Datalog[:, 4], Datalog[:, 1], "force", "pos_err", 2, 1)
    plot_uy(Datalog[:, 10], Datalog[:, 7], "interateout", "pid_put", 2, 1)
    u=chirp_cmd
    y=Datalog[:, 5]
    plot_uy(y, u, "ctrl_out", "chirp_cmd")#
    u=u-np.mean(u)
    y=y-np.mean(y)
    fw,transfer_fw=sys_id_correlate_resolution(gl_servo_freq, u, y, f_chirp_start/0.8, f_chirp_end/1.2,resolution=1.5)
    # fw,transfer_fw=sys_id_fft(gl_servo_freq, u, y, f_chirp_start/0.8, f_chirp_end/1.2)
    s_d = 1j * 2 * np.pi * fw
    z=np.e**(s_d*gl_T)
    s_z=1-z**(-1)#一阶后向差分,此处不除以T是因为Ki和Kd参数已经除过了
    c_fw=Kp+Ki*1/s_z+Kd*s_z
    close_loop_fw=np.ones_like(fw)-transfer_fw
    cp=close_loop_fw/(np.ones_like(fw)-close_loop_fw)
    pos_fw=cp/c_fw
    acc_fw=pos_fw*s_d*s_d
    bode_plot(fw,c_fw,"c_fw")
    bode_plot(fw,cp,"cp")
    bode_plot(fw,close_loop_fw,"close_loop_fw")
    bode_plot(fw,transfer_fw,"sensitivity")
    bode_plot(fw,acc_fw,"Acc")
    bode_plot(fw,pos_fw,"Pos")

    # 存储结果
    if enable_save_datalog:
        now_time = time.strftime('%Y_%m_%d_%H_%M_%S', time.localtime())
        path = os.getcwd()
        file_path = path + "\\Datalog\\" + now_time + '.csv'
        df1 = pd.DataFrame(data=Datalog,
                           columns=['sample', 'pos_err', 'set_pos', 'enc', 'force', 'ctrl_out', 'chirp_cmd',
                                    'pid_out', 'vel_fbk', 'set_jerk'])
        df1.to_csv(file_path, index=False)
def chirp_signal_plant_after_controller_close_loop_a_b(enable_save_datalog=False):
    plt.close("all")
    # 初始化扫频参数
    start_freq = 10
    stop_freq = 4000
    chirp_force=40
    T = 2#时间s
    f_chirp_start=start_freq*0.8
    f_chirp_end=stop_freq*1.2
    duration=(f_chirp_end-f_chirp_start)/(stop_freq-start_freq)*T
    t_list =np.arange(0,duration,1/gl_servo_freq)#step =1/servo_freq_current
    total_samples=len(t_list)
    chirp_cmd = chirp_force*np.sin(2*np.pi*((f_chirp_end-f_chirp_start)/duration*t_list**2/2+f_chirp_start*t_list))
    # 初始化控制器
    Kp, Ki, Kd = pid_design(method="hans_buttler", mass=mc.plant.mass, dampling=0.707, ratio=0.5,
                            alpha=mc.plant.c / mc.plant.mass, bandwidth_Hz=50)
    interateout = -20  #
    enc_fbk = 0  # 当前起始位置
    vel_fbk = 0  # 当前起始位置
    err = np.zeros(2)
    Datalog = np.zeros((total_samples, 10))
    profile=np.zeros((total_samples, 1))
    ffc_out=0
    # F_disturbance=np.sin(2*np.pi*15*t_list)#初始化外界干扰力
    F_disturbance=-1*interateout#初始化外界干扰力
    # F_disturbance=5*np.sin(2*np.pi*15*t_list)#初始化外界干扰力
    # F_disturbance=np.ones(total_samples)*25
    # snr = 10  # 信噪比
    # F_disturbance = wgn(F_disturbance, snr)  # 给构造电流采集信号加白噪声
    # begin runmng
    for i in range(total_samples):
        err[1] = profile[i] - enc_fbk
        Datalog[i, 1] = err[1]
        interateout += Ki * err[1]
        # 积分限幅
        interateout = 1950 if interateout > 1950 else interateout
        interateout = -1950 if interateout < -1950 else interateout

        pid_put = Kp * err[1] + interateout + Kd * (err[1] - err[0])
        # 控制器输出限幅
        pid_put = 2000 if pid_put > 2000 else pid_put
        pid_put = -2000 if pid_put < -2000 else pid_put
        err[0] = err[1]
        ctrl_out=pid_put+ffc_out+chirp_cmd[i]
        # 仿真实物运动
        # F_disturbance=10**5* enc_fbk+25 #模拟SSZ向上磁弹力
        force=ctrl_out+F_disturbance
        tout, y, x = ctrl.forced_response(gl_plant, T=[0, gl_T], U=[force, force],
                                          X0=[enc_fbk, vel_fbk])  # control库0.8.3
        vel_fbk = x[-1][-1]
        enc_fbk = y[-1]

        Datalog[i, 0] = i
        Datalog[i, 2] = profile[i]
        Datalog[i, 3] = enc_fbk
        Datalog[i, 4] = force
        Datalog[i, 5] = ctrl_out
        Datalog[i, 6] = chirp_cmd[i]
        Datalog[i, 7] = pid_put
        Datalog[i, 8] = vel_fbk
        Datalog[i, 9] = F_disturbance
    # 绘制结果
    plot_uy(Datalog[:, 9],chirp_cmd , "F_disturbance", "chirp_cmd", 2, 1)
    plot_uy(Datalog[:, 2], Datalog[:, 3], "set_pos", "enc", 2, 1)
    plot_uy(Datalog[:, 4], Datalog[:, 1], "force", "pos_err", 2, 1)
    plot_uy(Datalog[:, 7], Datalog[:, 5] , "pid_put", "chirp_cmd", 2, 1)
    u=Datalog[:, 5]+Datalog[:, 9]
    y=Datalog[:, 6]-Datalog[:, 5]
    # u=u-np.mean(u)
    # y=y-np.mean(y)
    plot_uy(u, y, "ctrl_out+Disturb", "chirp_cmd-ctrl_out")#所求为CP的传函

    # fw,transfer_fw=sys_id_correlate(gl_servo_freq, u, y, f_chirp_start/0.8, f_chirp_end/1.2)
    fw,transfer_fw=sys_id_correlate_resolution(gl_servo_freq, u, y, f_chirp_start/0.8, f_chirp_end/1.2,resolution=1.5)
    # fw,transfer_fw=sys_id_fft(gl_servo_freq, u, y, f_chirp_start/0.8, f_chirp_end/1.2)
    s_d = 1j * 2 * np.pi * fw
    z=np.e**(s_d*gl_T)
    s_z=1-z**(-1)#一阶后向差分,此处不除以T是因为Ki和Kd参数已经除过了
    c_fw=Kp+Ki*1/s_z+Kd*s_z
    # cp=(1-transfer_fw)/transfer_fw
    cp=transfer_fw
    pos_fw=cp/c_fw
    # close_loop_fw=1-transfer_fw
    close_loop_fw=cp/(1+cp)
    sensitivity=1/(1+cp)
    acc_fw=pos_fw*s_d*s_d
    bode_plot(fw,c_fw,"c_fw")
    bode_plot(fw,transfer_fw,"cp")
    bode_plot(fw,sensitivity,"sensitivity")
    bode_plot(fw,close_loop_fw,"close_loop_fw")
    bode_plot(fw,acc_fw,"Acc")
    bode_plot(fw,pos_fw,"Pos")

    # 存储结果
    if enable_save_datalog:
        now_time = time.strftime('%Y_%m_%d_%H_%M_%S', time.localtime())
        path = os.getcwd()
        file_path = path + "\\Datalog\\" + now_time + '.csv'
        df1 = pd.DataFrame(data=Datalog,
                           columns=['sample', 'pos_err', 'set_pos', 'enc', 'force', 'ctrl_out', 'chirp_cmd',
                                    'pid_out', 'vel_fbk', 'set_jerk'])
        df1.to_csv(file_path, index=False)
def chirp_signal_plant_after_controller_close_loop_c_d(enable_save_datalog=False):
    plt.close("all")
    # 初始化扫频参数
    start_freq = 30
    stop_freq = 4000
    chirp_force=40
    T = 1#时间s
    f_chirp_start=start_freq*0.8
    f_chirp_end=stop_freq*1.2
    duration=(f_chirp_end-f_chirp_start)/(stop_freq-start_freq)*T
    t_list =np.arange(0,duration,1/gl_servo_freq)#step =1/servo_freq_current
    total_samples=len(t_list)
    chirp_cmd = chirp_force*np.sin(2*np.pi*((f_chirp_end-f_chirp_start)/duration*t_list**2/2+f_chirp_start*t_list))
    # 初始化控制器
    Kp, Ki, Kd = pid_design(method="hans_buttler", mass=mc.plant.mass, dampling=0.707, ratio=0.5,
                            alpha=mc.plant.c / mc.plant.mass, bandwidth_Hz=50)
    interateout = -20  #
    enc_fbk = 0  # 当前起始位置
    vel_fbk = 0  # 当前起始位置
    err = np.zeros(2)
    Datalog = np.zeros((total_samples, 10))
    profile=np.zeros((total_samples, 1))
    ffc_out=0
    F_disturbance=-1*interateout#初始化外界干扰力
    # F_disturbance=5*np.sin(2*np.pi*15*t_list)#初始化外界干扰力
    # begin runmng
    for i in range(total_samples):
        err[1] = profile[i] - enc_fbk
        Datalog[i, 1] = err[1]
        interateout += Ki * err[1]
        # 积分限幅
        interateout = 1950 if interateout > 1950 else interateout
        interateout = -1950 if interateout < -1950 else interateout

        pid_put = Kp * err[1] + interateout + Kd * (err[1] - err[0])
        # 控制器输出限幅
        pid_put = 2000 if pid_put > 2000 else pid_put
        pid_put = -2000 if pid_put < -2000 else pid_put
        err[0] = err[1]
        ctrl_out=pid_put+ffc_out+chirp_cmd[i]
        # 仿真实物运动
        F_disturbance = -10 ** 3 * enc_fbk + 20  # 模拟SSZ向上磁弹力
        force=ctrl_out+F_disturbance
        tout, y, x = ctrl.forced_response(gl_plant, T=[0, gl_T], U=[force, force],
                                          X0=[enc_fbk, vel_fbk])  # control库0.8.3
        vel_fbk = x[-1][-1]
        enc_fbk = y[-1]

        Datalog[i, 0] = i
        Datalog[i, 2] = profile[i]
        Datalog[i, 3] = enc_fbk
        Datalog[i, 4] = force
        Datalog[i, 5] = ctrl_out
        Datalog[i, 6] = chirp_cmd[i]
        Datalog[i, 7] = pid_put
        Datalog[i, 8] = vel_fbk
        Datalog[i, 9] = F_disturbance
    # 绘制结果
    plot_uy(Datalog[:, 2], Datalog[:, 3], "set_pos", "enc", 2, 1)
    plot_uy(Datalog[:, 4], Datalog[:, 1], "force", "pos_err", 2, 1)
    u=chirp_cmd+Datalog[:, 9]
    y=Datalog[:, 5]+Datalog[:, 9]
    plot_uy(u, y, "chirp_cmd+disturb", "ctrl_out+disturb")#所求为CP的传函
    u=u-np.mean(u)
    y=u-np.mean(y)
    fw,transfer_fw=sys_id_correlate_resolution(gl_servo_freq, u, y, f_chirp_start/0.8, f_chirp_end/1.2,resolution=1.5)
    # fw,transfer_fw=sys_id_fft(gl_servo_freq, u, y, f_chirp_start/0.8, f_chirp_end/1.2)
    s_d = 1j * 2 * np.pi * fw
    z=np.e**(s_d*gl_T)
    s_z=1-z**(-1)#一阶后向差分,此处不除以T是因为Ki和Kd参数已经除过了
    c_fw=Kp+Ki*1/s_z+Kd*s_z
    close_loop_fw=np.ones_like(fw)-transfer_fw
    cp=close_loop_fw/(np.ones_like(fw)-close_loop_fw)
    pos_fw=cp/c_fw
    acc_fw=pos_fw*s_d*s_d
    bode_plot(fw,c_fw,"c_fw")
    bode_plot(fw,transfer_fw,"sensitivity")
    bode_plot(fw,close_loop_fw,"close_loop_fw")
    bode_plot(fw,acc_fw,"Acc")
    bode_plot(fw,pos_fw,"Pos")

    # 存储结果
    if enable_save_datalog:
        now_time = time.strftime('%Y_%m_%d_%H_%M_%S', time.localtime())
        path = os.getcwd()
        file_path = path + "\\Datalog\\" + now_time + '.csv'
        df1 = pd.DataFrame(data=Datalog,
                           columns=['sample', 'pos_err', 'set_pos', 'enc', 'force', 'ctrl_out', 'chirp_cmd',
                                    'pid_out', 'vel_fbk', 'set_jerk'])
        df1.to_csv(file_path, index=False)

def chirp_signal_plant_after_controller_close_loop_e_f(enable_save_datalog=False):
    plt.close("all")
    # 初始化扫频参数
    start_freq = 30
    stop_freq = 4000
    chirp_force=40
    T = 1#时间s
    f_chirp_start=start_freq*0.8
    f_chirp_end=stop_freq*1.2
    duration=(f_chirp_end-f_chirp_start)/(stop_freq-start_freq)*T
    t_list =np.arange(0,duration,1/gl_servo_freq)#step =1/servo_freq_current
    total_samples=len(t_list)
    chirp_cmd = chirp_force*np.sin(2*np.pi*((f_chirp_end-f_chirp_start)/duration*t_list**2/2+f_chirp_start*t_list))
    # 初始化控制器
    Kp, Ki, Kd = pid_design(method="hans_buttler", mass=mc.plant.mass, dampling=0.707, ratio=0.5,
                            alpha=mc.plant.c / mc.plant.mass, bandwidth_Hz=50)
    interateout = 0  #
    enc_fbk = 0  # 当前起始位置
    vel_fbk = 0  # 当前起始位置
    err = np.zeros(2)
    Datalog = np.zeros((total_samples, 10))
    profile=np.zeros((total_samples, 1))
    ffc_out=0
    F_disturbance=0#初始化外界干扰力
    F_disturbance=5*np.sin(2*np.pi*15*t_list)#初始化外界干扰力

    snr = 10  # 信噪比
    F_disturbance = wgn(F_disturbance, snr)  # 给构造电流采集信号加白噪声
    # begin runmng
    for i in range(total_samples):
        err[1] = profile[i] - enc_fbk
        Datalog[i, 1] = err[1]
        interateout += Ki * err[1]
        # 积分限幅
        interateout = 1950 if interateout > 1950 else interateout
        interateout = -1950 if interateout < -1950 else interateout

        pid_put = Kp * err[1] + interateout + Kd * (err[1] - err[0])
        # 控制器输出限幅
        pid_put = 2000 if pid_put > 2000 else pid_put
        pid_put = -2000 if pid_put < -2000 else pid_put
        err[0] = err[1]
        ctrl_out=pid_put+ffc_out+chirp_cmd[i]
        # 仿真实物运动
        # F_disturbance=25  #模拟SSZ向上磁弹力
        force=ctrl_out+F_disturbance[i]
        tout, y, x = ctrl.forced_response(gl_plant, T=[0, gl_T], U=[force, force],
                                          X0=[enc_fbk, vel_fbk])  # control库0.8.3
        vel_fbk = x[-1][-1]
        enc_fbk = y[-1]

        Datalog[i, 0] = i
        Datalog[i, 2] = profile[i]
        Datalog[i, 3] = enc_fbk
        Datalog[i, 4] = force
        Datalog[i, 5] = ctrl_out
        Datalog[i, 6] = chirp_cmd[i]
        Datalog[i, 7] = pid_put
        Datalog[i, 8] = vel_fbk
        Datalog[i, 9] = F_disturbance[i]
    # 绘制结果
    plot_uy(Datalog[:, 2], Datalog[:, 3], "set_pos", "enc", 2, 1)
    plot_uy(Datalog[:, 4], Datalog[:, 1], "force", "pos_err", 2, 1)
    u=chirp_cmd+Datalog[:, 9]
    y=-Datalog[:, 7]
    plot_uy(u, y, "chirp_cmd+disturb", "-pid_put")#所求为CP的传函
    plot_uy(Datalog[:, 9], Datalog[:, 5], "disturb", "ctrl_out")#所求为CP的传函
    # u=u-np.mean(u)
    # y=u-np.mean(y)
    fw,transfer_fw=sys_id_correlate_resolution(gl_servo_freq, u, y, f_chirp_start/0.8, f_chirp_end/1.2,resolution=1.5)
    # fw,transfer_fw=sys_id_fft(gl_servo_freq, u, y, f_chirp_start/0.8, f_chirp_end/1.2)
    s_d = 1j * 2 * np.pi * fw
    z=np.e**(s_d*gl_T)
    s_z=1-z**(-1)#一阶后向差分,此处不除以T是因为Ki和Kd参数已经除过了
    c_fw=Kp+Ki*1/s_z+Kd*s_z
    close_loop_fw=transfer_fw
    cp=close_loop_fw/(np.ones_like(fw)-close_loop_fw)
    pos_fw=cp/c_fw
    acc_fw=pos_fw*s_d*s_d
    bode_plot(fw,c_fw,"c_fw")
    bode_plot(fw,cp,"cp")
    bode_plot(fw,close_loop_fw,"close_loop_fw")
    bode_plot(fw,acc_fw,"Acc")
    bode_plot(fw,pos_fw,"Pos")

    # 存储结果
    if enable_save_datalog:
        now_time = time.strftime('%Y_%m_%d_%H_%M_%S', time.localtime())
        path = os.getcwd()
        file_path = path + "\\Datalog\\" + now_time + '.csv'
        df1 = pd.DataFrame(data=Datalog,
                           columns=['sample', 'pos_err', 'set_pos', 'enc', 'force', 'ctrl_out', 'chirp_cmd',
                                    'pid_out', 'vel_fbk', 'set_jerk'])
        df1.to_csv(file_path, index=False)
if __name__ == '__main__':
    dB=-22.2879
    print(1/10**(dB/20))
    set_machine_constant()
    k=1640#10**5
    f=np.sqrt(k/13)/2/np.pi
    print("k,f=%s"%[k,f])
    # run_closed_loop()
    # chirp_signal_plant_open_loop()
    # chirp_signal_plant_open_loop_acc()

    # chirp_signal_plant_after_controller_close_loop_sigin_pidout(enable_save_datalog=False)#1
    chirp_signal_plant_after_controller_close_loop_sigin_ctrlout(enable_save_datalog=False)#2
    # chirp_signal_plant_after_controller_close_loop_a_b(enable_save_datalog=False)#3
    # chirp_signal_plant_after_controller_close_loop_c_d(enable_save_datalog=False)#4
    # chirp_signal_plant_after_controller_close_loop_e_f(enable_save_datalog=False)#5
    plt.show()
    plt.close("all")
    # disturb_pa_noise = np.random.random(len(t_list)) / 100
    # quick_plot(t=t_list, yt=disturb_pa_noise, title="disturb_pa_noise")



