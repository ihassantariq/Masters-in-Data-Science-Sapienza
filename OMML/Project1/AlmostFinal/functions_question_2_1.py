import numpy as np
import time
from scipy.optimize import minimize
from sklearn.model_selection import train_test_split
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import codecs,json

def data_set():
    df = pd.read_csv("DATA.csv")
    X = np.array(df[['x1','x2']])
    y = np.array(df[['y']])
    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.75)
    return X_train.T, y_train.T, X_test.T, y_test.T

def tanh(x,sigma):
    return (1-np.exp(-2*x*sigma))/(1+np.exp(-2*x*sigma))

def Xpredict(x, w, b, v,sigma):
    return np.dot(v, tanh(np.dot(w, x) + b, sigma))

def loss_function(v, w,b, X, y,rho,sigma):
    return (1/2)*np.mean(np.square(Xpredict(X,w,b,v,sigma) - y) ) + rho*np.linalg.norm(v,ord=2)**2

def Xfit(X, y,w, b, v, rho, sigma):
    t = time.time()
    res = minimize(loss_function, x0=v, args=(w, b, X, y, rho, sigma),method='BFGS')
    tr_Time = time.time()-t
    return w, b, res.x, res.nfev, res.njev, tr_Time

def MSE(w,b,v,x,y,sigma):
    return 0.5*np.mean(np.square(Xpredict(x,w,b,v,sigma) - y))
    
def XLearn(X_tr, y_tr, X_ts, y_ts, neurons, rho, sigma, input_dim, output_dim, gens):
    error=[]
    terror=[]
    test_error = 1000
    for times in range(gens):
        np.random.seed(1870543+times)
        w = np.random.randn(neurons, input_dim)
        b = np.random.randn(neurons, 1)
        v = np.random.randn(output_dim, neurons) 

        w_opt, b_opt, v_opt, FE_opt, GE_opt,tr_Time = Xfit(X_tr,y_tr,w,b,v,rho,sigma)
        error.append(MSE(w_opt,b_opt,v_opt,X_ts,y_ts,sigma))
        terror.append(MSE(w_opt,b_opt,v_opt,X_tr,y_tr,sigma))
        if error[-1] < test_error:
            test_error = error[-1]
            train_error = terror[-1]
            train_time = tr_Time
            w = w_opt
            v = v_opt
            b = b_opt
            FE = FE_opt
            GE = GE_opt
    return test_error, train_error, train_time, w, v, b, FE, GE, error, terror

def writejson(file, w, b, v, neurons, sigma, rho):
    values = {'w':w.tolist(),'b':b.tolist(),'v':v.tolist(),'neurons':neurons,'sigma':sigma,'rho':rho}
    json.dump(values, codecs.open(file, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=False, indent=4)

def plotting(file,w,b,v,sigma):
    fig = plt.figure(figsize=(12,8))
    ax = plt.axes(projection='3d')
    x = np.linspace(-2, 2, 1000) 
    y = np.linspace(-1, 1, 1000)
    X, Y = np.meshgrid(x, y)
    zs = Xpredict(np.array([np.ravel(X),np.ravel(Y)]),w,b,v,sigma)
    Z= zs.reshape(X.shape)
    ax.plot_surface(Y, X, Z, cmap='viridis', edgecolor='none')
    ax.set_xlabel('y')
    ax.set_ylabel('x')
    ax.set_zlabel('z')
    plt.savefig(file,dpi=600)
    plt.show()
    
def writetxt(file, neurons, sigma, rho, method, t_error, min_error, tr_time, FE, GE):
    print('\nSaving results in output_question_2_1.txt')
    output = open(file,"w")
    output.write("This is homework 1: question 1_1")
    output.write("\nNumber of neurons: " + "%i" % neurons)
    output.write("\nsigma parameter: " + "%.1f" % sigma)
    output.write("\nrho parameter: " + "%f" % rho)
    output.write("\nMethod:"+method)
    output.write("\nFunction evaluations: " + "%i" % FE)
    output.write("\nGradient evaluations: " + "%i" % GE)
    output.write("\nTraining time: " + "%f" % tr_time + " sec")
    output.write("\nTraining error: " + "%f" % t_error)
    output.write("\nTest error: " + "%f" % min_error)
    output.close()

def ovunfit(gens,error,terror):
    plt.plot([i for i in range(gens)],error,label='test_error')
    plt.plot([i for i in range(gens)],terror,label='train_error')
    plt.legend()
    plt.xlabel('Number of generations')
    plt.ylabel('Error Values')
    plt.title('Error Dynamics')
    plt.grid(True)
    plt.savefig('error_question_2_1.png',dpi=600)
    plt.show()