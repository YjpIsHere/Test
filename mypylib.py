import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft,ifft,fftshift
import scipy
import scipy.signal as signal
import pandas as pd
import os

def correlate_full(a,b):
    #等价于np.correlate(a,b,"full")
    num = len(a)
    left_pad = round((num-1)/2+0.1)
    right_pad = num-1-left_pad
    a1 = np.pad(a,(left_pad,right_pad),"constant",constant_values=0)
    b1 = np.pad(np.flip(b),(left_pad,right_pad),"constant",constant_values=0)
    ahat = fft(a1)
    bhat = fft(b1)
    chat = ahat*bhat
    c1 = np.real(ifft(chat))
    c1_fftshift = fftshift(c1)
    if num%2==1:
        c1_fftshift = np.roll(c1_fftshift,1)
    return c1_fftshift


def sys_id_correlate(fs_Hz,u,y,start_freq,stop_freq):
    #resolution is defined by fs_Hz/M
    num = len(y)
    if fs_Hz == 1000000:
        M = 4096
    else:
        M = 1024
    resolution = 1.5#fs_Hz/(num*2)
    Ruu = np.zeros(2*M)
    Ryu = np.zeros(2*M)
    Ruy = np.zeros(2*M)
    ruu = correlate_full(u,u)
    ruy = correlate_full(u,y)
    ryu = correlate_full(y,u)

    if M>num:
        Ruu[:num] = np.flip(ruu[:num])
        Ruy[:num] = np.flip(ruy[:num])
        Ryu[:num] = np.flip(ryu[:num])
    else:
        Ruu[:M] = np.flip(ruu[num-M:num])
        Ruy[:M] = np.flip(ruy[num-M:num])
        Ryu[:M] = np.flip(ryu[num-M:num])

    Ruu = Ruu/num
    Ruy = Ruy/num
    Ryu = Ryu/num
    window = 0.5*(1+np.cos(2*np.pi*np.arange(M)/M))
    Ruu[0:M] = Ruu[0:M]*window
    Ruy[0:M] = Ruy[0:M]*window
    Ryu[0:M] = Ryu[0:M]*window

    Ruu[M+1:] = np.flip(Ruu[1:M])
    Ruy[M+1:] = np.flip(Ryu[1:M])

    start_point = int(start_freq/resolution)
    stop_point = int(stop_freq/resolution)+1
    f = np.arange(start_point,stop_point)*resolution
    fRuu = fft(Ruu)[start_point:stop_point]
    fRyu = fft(Ruy)[start_point:stop_point]
    return f,fRyu/fRuu

def sys_id_correlate_resolution(fs_Hz,u,y,start_freq,stop_freq,resolution):
    #resolution is defined by user
    num = len(y)
    M=int(fs_Hz/resolution/2)
    Ruu = np.zeros(2*M)
    Ryu = np.zeros(2*M)
    Ruy = np.zeros(2*M)
    ruu = correlate_full(u,u)
    ruy = correlate_full(u,y)
    ryu = correlate_full(y,u)

    if M>num:
        Ruu[:num] = np.flip(ruu[:num])
        Ruy[:num] = np.flip(ruy[:num])
        Ryu[:num] = np.flip(ryu[:num])
    else:
        Ruu[:M] = np.flip(ruu[num-M:num])
        Ruy[:M] = np.flip(ruy[num-M:num])
        Ryu[:M] = np.flip(ryu[num-M:num])

    Ruu = Ruu/num
    Ruy = Ruy/num
    Ryu = Ryu/num
    window = 0.5*(1+np.cos(2*np.pi*np.arange(M)/M))
    Ruu[0:M] = Ruu[0:M]*window
    Ruy[0:M] = Ruy[0:M]*window
    Ryu[0:M] = Ryu[0:M]*window

    Ruu[M+1:] = np.flip(Ruu[1:M])
    Ruy[M+1:] = np.flip(Ryu[1:M])

    start_point = int(start_freq/resolution)
    stop_point = int(stop_freq/resolution)+1
    f = np.arange(start_point,stop_point)*resolution
    fRuu = fft(Ruu)[start_point:stop_point]
    fRyu = fft(Ruy)[start_point:stop_point]
    return f,fRyu/fRuu
def sys_id_fft(fs_Hz,u,y,start_freq,stop_freq):
    M = len(y)
    resolution = fs_Hz/M
    start_point=int(start_freq/resolution)
    stop_point=int(stop_freq/resolution)
    f = np.arange(start_point,stop_point)*resolution
    fRuu = fft(u)[start_point:stop_point]
    fRyu = fft(y)[start_point:stop_point]
    return f,fRyu/fRuu


def bode_plot(f_Hz,transfer_function_fw,title):
    gain = 20 * np.log10(np.abs(transfer_function_fw))
    phase = np.angle(transfer_function_fw, deg=True)
    fig, axes = plt.subplots(2, 1)
    fig.suptitle(title)
    axes[0].set_xlabel("frequency(Hz)")
    axes[0].set_ylabel("Gain(dB)")
    axes[0].semilogx(f_Hz, gain)
    axes[0].grid()

    axes[1].set_xlabel("frequency(Hz)")
    axes[1].set_ylabel("Phase(deg)")
    axes[1].semilogx(f_Hz, phase)
    axes[1].grid()
    plt.tight_layout()
    plt.savefig("bode.png")

def plot_uy(u,y,u_title,y_title,row=None,col=None):
    if row == None:
        plt.figure()
        plt.plot(u,label=u_title)
        plt.plot(y,label=y_title)
        plt.title('<'+u_title+'>'+' and '+'<'+y_title+'>')
        plt.legend(loc='best')
        plt.tight_layout()
        plt.grid()
    else:
        fig, axes = plt.subplots(row, col)
        plt.title('<' + u_title + '>' + ' and ' + '<' + y_title + '>')
        axes[0].plot(u,color='r',label=u_title)
        axes[0].grid()
        plt.grid()

        axes[1].plot(y,color='g',label=y_title)
        axes[1].grid()
        plt.grid()
    plt.savefig(u_title+" and "+y_title+".png")

def quick_plot(t=None,yt=None,title=""):
    plt.figure()
    if t is None:
        plt.plot(yt,label=title)
    else:
        plt.plot(t,yt,color='r',label=title)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.grid()

def DFT_single_freq(sample_freq,freq,signal):
    '''
    输入：
        sample_freq：采样频率
        signal：指输入信号，我希望它是array格式
        freq：需要计算单一频率DFT的频率，单位Hz
    输出：
        单一频率DFT计算出的结果,是个复数
    '''
    out_put=0+0j
    T=1/sample_freq
    for i in range(len(signal)):
        w=-2*np.pi*freq*i*T
        out_put+=signal[i]*np.e**(1j*w)
    return out_put

def wgn(x, snr):
    len_x = len(x)
    Ps = np.sum(np.power(x, 2)) / len_x
    Pn = Ps / (np.power(10, snr / 10))
    noise = np.random.randn(len_x) * np.sqrt(Pn)
    return x + noise
