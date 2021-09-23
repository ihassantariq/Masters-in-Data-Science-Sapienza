import numpy as np
import time
from scipy.optimize import minimize
from sklearn.model_selection import train_test_split
import pandas as pd
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import codecs,json
'''
This function is needed to split the data to train and test sets
'''
def data_set():
    df = pd.read_csv("DATA.csv")
    X = np.array(df[['x1','x2']])
    y = np.array(df[['y']])
    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.85)
    return X_train.T, y_train.T, X_test.T, y_test.T

#
def tanh(x,sigma):
    return (1-np.exp(-2*x*sigma))/(1+np.exp(-2*x*sigma))

def wrap(w, b, v):
    return np.concatenate([w.reshape(-1),b.reshape(-1),v.reshape(-1)])

def unwrap(inD, neu, outD, wraped):
    w = wraped[0:inD*neu].reshape(neu,inD)
    b = wraped[inD*neu:inD*neu+outD*neu].reshape(neu,outD)
    v = wraped[inD*neu+outD*neu:inD*neu+outD*neu + neu*outD].reshape(outD,neu)
    return w,b,v

def predict(x, w, b, v,sigma):
    return np.dot(v, tanh(np.dot(w, x) + b, sigma))

def loss_function(params,X,y,rho,sigma,input_dim,neurons):
    w,b,v = unwrap(input_dim,neurons,1,params)
    return (1/2)*np.mean( np.square(predict(X, w, b, v, sigma) - y) ) + rho*np.linalg.norm(params, ord=2)**2

def fit(X, y,neurons, w, b, v, rho, sigma, input_dim):
    t = time.time()
    res = minimize(loss_function, x0=wrap(w, b, v), args=(X, y, rho, sigma, input_dim, neurons),method='BFGS', options={'disp':False})
    tr_Time = time.time()-t
    w, b, v = unwrap(2, neurons, 1, res.x)
    FE = res.nfev
    GE = res.njev
    return w, b, v, FE, GE,tr_Time

def MSE(w,b,v,x,y,sigma):
    return 0.5*np.mean(np.square(predict(x,w,b,v,sigma) - y))

def K_fold(n_folds, X_tr,y_tr,neurons,w,b,v,rho,sigma, input_dim):
                for i in range(n_folds):
                    X_train, X_val, y_train, y_val = train_test_split(X_tr.T, y_tr.T, train_size=0.8, shuffle = True)
                    w, b, v, FE, GE,tr_Time = fit(X_train.T,y_train.T,neurons,w,b,v,rho,sigma,input_dim)
                    #print('CV-error', MSE(w,b,v,X_val.T,y_val.T,sigma))
                return w, b, v, FE, GE, tr_Time

def gridsearch(X_tr, y_tr, X_ts, y_ts, input_dim,output_dim,N,S,R):
    error = []
    terror = []
    test_error = 100
    for sigma in S:
        for rho in R:
            for neurons in N: 
                w = np.random.randn(neurons, input_dim)
                b = np.random.randn(neurons, 1)
                v = np.random.randn(output_dim, neurons)
                
                w_opt,b_opt,v_opt,FE_opt,GE_opt,tr_Time = K_fold(5, X_tr,y_tr,neurons,w,b,v,rho,sigma, input_dim)
                #w_opt,b_opt,v_opt,FE_opt,GE_opt,tr_Time = fit(X_tr,y_tr,neurons,w,b,v,rho,sigma,input_dim)
                error.append(MSE(w_opt,b_opt,v_opt,X_ts,y_ts,sigma))
                terror.append(MSE(w_opt,b_opt,v_opt,X_tr,y_tr,sigma))
                if error[-1] < test_error:
                    test_error = error[-1]
                    train_error = terror[-1]
                    train_time = tr_Time
                    n = neurons
                    s = sigma
                    r = rho
                    w = w_opt
                    v = v_opt
                    b = b_opt
                    FE = FE_opt
                    GE = GE_opt
    return test_error, train_error, train_time, n, s, r, w, v, b,FE,GE, error, terror

def writetxt(file, neurons, sigma, rho, method, t_error, min_error, tr_time, FE, GE):
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
    output.write("\nTest error: " + "%f" % min_error)
    output.close()
    
def writejson(file, w_opt, b_opt, v_opt, neurons, sigma, rho):
    values = {'w':w_opt.tolist(),'b':b_opt.tolist(),'v':v_opt.tolist(),'neurons':neurons,'sigma':sigma,'rho':rho}
    file_path = file 
    json.dump(values, codecs.open(file_path, 'w', encoding='utf-8'), separators=(',', ':'), sort_keys=False, indent=4)


def plotting(file,w,b,v,sigma):
    fig = plt.figure(figsize=(12,8))
    ax = plt.axes(projection='3d')
    x = np.linspace(-2, 2, 1000) 
    y = np.linspace(-1, 1, 1000)
    X, Y = np.meshgrid(x, y)
    zs = predict(np.array([np.ravel(X),np.ravel(Y)]),w,b,v,sigma)
    Z= zs.reshape(X.shape)
    ax.plot_surface(Y, X, Z, cmap='viridis', edgecolor='none')
    ax.set_xlabel('y')
    ax.set_ylabel('x')
    ax.set_zlabel('z')
    plt.savefig(file,dpi=600)
    plt.show()
    

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
    plt.savefig('error_question_1_1.png',dpi=600)
    plt.show()