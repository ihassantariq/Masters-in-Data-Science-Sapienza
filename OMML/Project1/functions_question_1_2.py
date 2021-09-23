import time
import numpy as np
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
    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.85)
    return X_train.T, y_train.T, X_test.T, y_test.T

def gauss(x,c,sigma):
    return np.exp(-np.square((np.linalg.norm(x-c, ord=1, axis=1)/sigma)))

def wrap(c,v):
    return np.concatenate([c.reshape(-1),v.reshape(-1)])

def unwrap(inD, neu, outD, wraped):
    c = wraped[0:inD*neu].reshape((inD,neu))
    v = wraped[inD*neu:inD*neu + outD*neu].reshape((outD,neu))
    return c,v

def loss_function(params,X,y,rho,sigma,input_dim,neurons):
    c,v = unwrap(input_dim,neurons,1,params)
    return (1/2)*np.mean( np.square(RBFpredict(X,c,v,sigma) - y) ) + rho*np.linalg.norm(params, ord=2)**2

def RBFpredict(x, c, w, sigma):
    return np.dot(w, np.array([gauss(x.T, c[:,c_idx], sigma) for c_idx in range(c.shape[1])]).reshape(c.shape[1],x.shape[1]))

def RBFfit(X,y,neurons, c, v, rho, sigma, input_dim):
    t = time.time()
    res = minimize(loss_function, x0=wrap(c, v), args=(X, y, rho, sigma, input_dim, neurons),method='BFGS')
    training_Time = time.time()-t
    c_opt, v_opt = unwrap(2, neurons, 1, res.x)
    return c_opt, v_opt, res.fun, res.nit, res.nfev, res.njev, training_Time, res

def K_fold(n_folds, X_tr,y_tr,neurons,c,v,rho,sigma, input_dim):
    for i in range(n_folds):
        X_train, X_val, y_train, y_val = train_test_split(X_tr.T, y_tr.T, train_size=0.8, shuffle = True)
        c_opt,v_opt,func_opt,num_Iter,FE_opt,GE_opt,tr_Time,res = RBFfit(X_train.T,y_train.T,neurons,c,v,rho,sigma,input_dim)
        #print('CV-error', MSE(c,v,X_val.T,y_val.T,sigma))
    return c_opt, v_opt, func_opt, num_Iter, FE_opt, GE_opt, tr_Time, res

def gridsearch(X_tr, y_tr, X_ts, y_ts, input_dim,output_dim,N,S,R):
    error = []
    terror = []
    test_error = 100
    for sigma in S:
        for rho in R:
            for neurons in N:
                c = np.random.randn(2,neurons)
                v = np.random.randn(1,neurons)
                c_opt,v_opt,func_opt,num_Iter,FE_opt,GE_opt,tr_Time,res = K_fold(5, X_tr,y_tr,neurons,c,v,rho,sigma, input_dim)
                
                error.append(MSE(c_opt, v_opt, X_ts, y_ts, sigma))
                terror.append(MSE(c_opt,v_opt,X_tr,y_tr,sigma))
                
                if error[-1] < test_error:
                    test_error = error[-1]
                    train_error = terror[-1]
                    train_time = tr_Time
                    n = neurons
                    s = sigma
                    r = rho
                    c = c_opt
                    v = v_opt
                    FE = FE_opt
                    GE = GE_opt
    return test_error, train_error, train_time, n, s, r, c, v, FE, GE, error, terror
    
def MSE(c,v,x,y,sigma):
    return 0.5*np.mean((RBFpredict(x,c,v,sigma) - y) **2)

def plotting(file,c,v,sigma):
    fig = plt.figure(figsize=(12,8))
    ax = plt.axes(projection='3d')
    x = np.linspace(-2, 2, 1000) 
    y = np.linspace(-1, 1, 1000)
    X, Y = np.meshgrid(x, y)
    zs = RBFpredict(np.array([np.ravel(X),np.ravel(Y)]),c,v,sigma)
    Z= zs.reshape(X.shape)
    ax.plot_surface(Y, X, Z, cmap='viridis', edgecolor='none')
    ax.set_xlabel('y')
    ax.set_ylabel('x')
    ax.set_zlabel('z')
    plt.savefig(file,dpi=600)
    plt.show()

def writetxt(file, neurons, sigma, rho, method, t_error, test_error, tr_time, FE, GE):
    print('\nSaving results in output_question_1_1.txt')
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
    output.write("\nTest error: " + "%f" % test_error)
    output.close()
    
def writejson(file, c_opt, v_opt, neurons, sigma, rho):
    values = {'c':c_opt.tolist(),'v':v_opt.tolist(),'neurons':neurons,'sigma':sigma,'rho':rho}
    file_path = file 
    json.dump(values, codecs.open(file_path, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=False, indent=4)
    
def ovunfit(N,S,R,sigma,rho,error,terror):
    plt.plot(N,np.array(error).reshape(len(S),len(R),len(N))[S.index(sigma)][R.index(rho)],label='test_error')
    plt.plot(N,np.array(terror).reshape(len(S),len(R),len(N))[S.index(sigma)][R.index(rho)],label='train_error')
    plt.legend()
    plt.xlabel('Number of neurons')
    plt.ylabel('Error Values')
    plt.title('Error Dynamics')
    plt.text(N[-1]-1, error[-1]/2, r'$\rho=$'+str(rho))
    plt.text(N[-1]-1, terror[-1]/3, r'$\sigma=$'+str(sigma))
    plt.grid(True)
    plt.savefig('error_question_1_2.png',dpi=600)
    plt.show()