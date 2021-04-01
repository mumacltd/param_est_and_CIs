import numpy as np
import json
import unittest


coneThreshold = -2
coneOffset = 2
coneTimeConstant = 1
rodS2 = -0.1
rodAlpha = 5

parameters = np.array(
    [coneThreshold, coneOffset, coneTimeConstant, rodS2, rodAlpha])


gradientVectorReference = np.array(
    [-12.6998098, - 1.8778356, - 5.8712005,  17.6243213,   0.2213682])

time = np.array([0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75,
                 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6,
                 6.25, 6.5, 6.75, 7, 7.25, 7.5, 7.75, 8, 8.25, 8.5, 8.75, 9, 9.25,
                 9.5, 9.75, 10])

threshold = np.array([-0.460309841006592, -0.620877386262897,
                      -0.729923060588993, -0.894630651754359, -1.02666933827986, -1.11357313776059,
                      -1.22467356210846, -1.19842194849487, -1.34572163740579, -1.36800467390198,
                      -1.39102152095425, -1.49690208693398, -1.45249997847129, -1.49870799087709,
                      -1.5351248783815, -1.53613994152801, -1.58086202577915, -1.6487659231982,
                      -1.5884389355972, -1.69341314952358, -1.69970119019474, -1.83596425529381,
                      -1.85320861021247, -1.91939817342708, -1.99943092645119, -2.09344257568903,
                      -2.11389605374682, -2.21941294749368, -2.2618972051465, -2.31991456775266,
                      -2.31949884431587, -2.4638563047387, -2.49960126740859, -2.57274990082888,
                      -2.57800015816419, -2.63261915158418, -2.65289491063321, -2.7852209514604,
                      -2.77068149386663, -2.86272031406487])


def Heaviside(time, alpha):
    return np.round(1/(1 + np.exp(-20*(time - alpha))), 0)


def coneModel(parameters, time):
    Cone = parameters[0] + parameters[1] * np.exp(-time/parameters[2])
    return Cone


def rodModel(parameters, time):
    Rod = parameters[3] * (time - parameters[4]) * \
        Heaviside(time, parameters[4])
    return Rod


def fiveParameterModel(parameters, time):
    Cone = coneModel(parameters, time)
    Rod1 = rodModel(parameters, time)
    Val = Rod1 + Cone
    return(Val)


def residuals(parameters, time, threshold):
    Val = threshold - fiveParameterModel(parameters, time)
    return Val


def leastSquares(parameters, time, threshold):
    a = residuals(parameters, time, threshold)
    Val = np.dot(a, a)
    return Val


def coneThresholdPDE(parameters, time, threshold):
    root = time/time
    Val = -2*np.dot(root, residuals(parameters, time, threshold))
    return(Val)


def coneOffsetPDE(parameters, time, threshold):
    root = np.exp(-time/parameters[2])
    Val = -2*np.dot(root, residuals(parameters, time, threshold))
    return(Val)


def coneTimeConstantPDE(parameters, time, threshold):
    root = np.exp(-time/parameters[2])*parameters[1] * \
        time/(parameters[2]*parameters[2])
    Val = -2*np.dot(root, residuals(parameters, time, threshold))
    return(Val)


def rodS2PDE(parameters, time, threshold):
    root = (time-parameters[4])/(1 + np.exp(-20*(time-parameters[4])))
    Val = -2*np.dot(root, residuals(parameters, time, threshold))
    return(Val)


def rodAlphaPDE(parameters, time, threshold):
    k = 20
    expr6 = time - parameters[4]
    expr7 = parameters[3] * expr6
    expr10 = np.exp(-k * expr6)
    expr11 = 1 + expr10
    root = (parameters[3]/expr11 + expr7 * (expr10 * k)/np.power(expr11, 2))
    Val = 2*np.dot(root, residuals(parameters, time, threshold))
    return(Val)


def gradientVector(parameters, time, threshold):
    Val = np.array([coneThresholdPDE(parameters, time, threshold),
                    coneOffsetPDE(parameters, time, threshold),
                    coneTimeConstantPDE(parameters, time, threshold),
                    rodS2PDE(parameters, time, threshold),
                    rodAlphaPDE(parameters, time, threshold)])
    return Val


tmp = np.round(gradientVector(parameters, time, threshold) /
               gradientVectorReference, 5)

print(tmp)

initialParameters = parameters


def findParametersGradDescentPlain(parameters, time, threshold):
    tol = 1e-6
    eta = 1/np.array([1000, 1000, 1000, 10000, 10])
    LS0 = 1e7
    for ii in range(10001):
        gradient = gradientVector(parameters, time, threshold)
        parameters = parameters - eta*gradient
        LS1 = leastSquares(parameters, time, threshold)
        delta = np.abs(LS0 - LS1)
        LS0 = LS1
        if delta < tol:
            break
    Val = {'parameters': parameters, 'runs': ii,
           'noise': np.std(residuals(parameters, time, threshold))}
    return Val


# findParametersGradDescentPlain(initialParameters, time, threshold)
# {'parameters': array([-1.61953284,  1.43606618,  1.41053068, -0.22708588,  4.47991306]),
# 'runs': 3813,
# 'noise': 0.033457640158571667}


def findParametersGradDescentMomentum(parameters, time, threshold):
    tol = 1e-6
    eta = 1/np.array([1000, 1000, 1000, 10000, 10])
    LS0 = 1e7
    velocity = 0
    gamma = 0.9
    for ii in range(10001):
        gradient = gradientVector(parameters, time, threshold)
        velocity = gamma * velocity + eta*gradient
        parameters = parameters - velocity
        LS1 = leastSquares(parameters, time, threshold)
        delta = np.abs(LS0 - LS1)
        LS0 = LS1
        if delta < tol:
            break
    Val = {'parameters': parameters, 'runs': ii,
           'noise': np.std(residuals(parameters, time, threshold))}
    return Val


# findParametersGradDescentMomentum(initialParameters, time, threshold)
# {'parameters': array([-1.69372693,  1.45025906,  1.66935426, -0.22427115,  4.72741029]),
# 'runs': 819,
# 'noise': 0.02959918420056771}


def findParametersGradDescentNAG(parameters, time, threshold):
    tol = 1e-6
    eta = 1/np.array([1000, 1000, 1000, 10000, 10])
    LS0 = 1e7
    velocity = 0
    gamma = 0.9
    for ii in range(10001):
        gradient = gradientVector(
            parameters - velocity * gamma, time, threshold)
        velocity = gamma * velocity + eta*gradient
        parameters = parameters - velocity
        LS1 = leastSquares(parameters, time, threshold)
        delta = np.abs(LS0 - LS1)
        LS0 = LS1
        if delta < tol:
            break
    Val = {'parameters': parameters, 'runs': ii,
           'noise': np.std(residuals(parameters, time, threshold))}
    return Val


# findParametersGradDescentNAG(initialParameters, time, threshold)
# {'parameters': array([-1.69370684,  1.45026001,  1.66926942, -0.22427364,  4.72737398]), 
# 'runs': 1414, 
# 'noise': 0.029599570676824342}
