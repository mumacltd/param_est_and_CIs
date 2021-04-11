// var out = null;
// var parameters = [15.130680495, -13.579104620999999, 0.6162903851, 1.751947972, 3.0147202696];

// var rawData = [{x: 0.2629666666666667, y: "6.13"}, {x: 0.4753, y: "8.86"}, {x: 0.7072333333333334, y: "10.14"}, {x: 0.9555, y: "12.45"}, {x: 1.2155333333333334, y: "13.59"}, {x: 1.4863333333333333, y: "13.82"}, {x: 1.7687166666666667, y: "15.05"}, {x: 2.054666666666667, y: "14.15"}, {x: 2.3568, y: "15.62"}, {x: 2.6583, y: "14.73"},]
// parameters = [-1.47307,  1.38009,  0.58676, -0.21387,  2.92516 ]
 //obj = [
 //   [2.629667e-01,4.753000e-01,7.072333e-01,9.555000e-01,1.215533e+00,1.486333e+00,1.768500e+00,2.055267e+00,2.348800e+00,2.652500e+00,2.953783e+00,3.245683e+00,3.524033e+00,3.806617e+00,4.102017e+00,4.394400e+00,4.684317e+00,4.972533e+00],
 //   [-6.140000e-01,-8.870000e-01,-9.270000e-01,-1.247000e+00,-1.361000e+00,-1.383000e+00,-1.390000e+00,-1.461000e+00,-1.409000e+00,-1.442000e+00,-1.394000e+00,-1.529000e+00,-1.638000e+00,-1.768000e+00,-1.666000e+00,-1.816000e+00,-1.879000e+00,-1.839000e+00]]
//


function extractData(rawData){
    
    var time = rawData.map(a => a.x);
    var str = rawData.map(a=> a.y);
    var thrs = str.map(Number);

    obj = [time, thrs];
return obj;
}
//var obj = extractData(rawData);

//console.log(obj)


function heaviside(parameters, obj){
    
    var val = parameters[4];
    var time = obj[0];
    out = time.map(e => Math.round(1/(1+ Math.exp(-20*(e-val)))));
    
    return out;
}

// // console.log(heaviside(parameters, obj))

function coneModel(parameters, obj){
    ct = parameters[0];
    cc = parameters[1];
    tau = parameters[2];
    time =obj[0];
    out = time.map(e => ct + cc*Math.exp(-e/tau));
    return out;

}
// // console.log(coneModel(parameters, obj));

function rodModel(parameters, obj){
    time = obj[0];
    S2 = parameters[3];
    alpha = parameters[4]
    a = time.map(e => S2 * (e - alpha));
    b = heaviside(parameters,obj);
    out = a.map((e,i)=>e * b[i])
    
    return out;

}

function fiveParameterModel(parameters, obj){
    cone = coneModel(parameters, obj);
    rod = rodModel(parameters, obj);
    out = cone.map((e,i) => e + rod[i])

    return out;
}
// console.log(fiveParameterModel(parameters, obj))

function residuals(parameters, obj){
    thrs = obj[1];
    model = fiveParameterModel(parameters,obj)

    out = thrs.map((e,i) => e - model[i])
    return out
}
// console.log(residuals(parameters, obj));

function dot(a,b){
    out = a.map((e,i) => e *b[i]);
    out = out.reduce((m,n)=> m+n,0)

    return out; 
}
// console.log(residuals(parameters, obj))

function leastSquare(parameters, obj){
    a = residuals(parameters, obj);
    out = a.map((e,i) => e *a[i]);
    out = out.reduce((m,n)=> m+n,0)

    return out; 
}
// console.log(leastSquare(parameters, obj))

function coneThresholdPDE(parameters, obj){
    time = obj[0]
    root = time.map(e => e/e);
    resid = residuals(parameters,obj)

    out = -2*dot(root,resid)
    return out
}
// out =coneThresholdPDE(parameters, obj)
//console.log(coneThresholdPDE(parameters, obj))


function coneThresholdPDEVec(parameters, obj){
    time = obj[0]
    root = time.map(e => e/e);

    return root
}
// out = coneThresholdPDEVec(parameters, obj)
// console.log(coneThresholdPDEVec(parameters, obj))

function coneCoeffPDE(parameters, obj){
    time = obj[0]
    CC = parameters[2];
    root = time.map(t => Math.exp(-t/CC));
    resid = residuals(parameters,obj)

    out = -2*dot(root,resid)
    return out
}



function coneCoeffPDEVec(parameters, obj){
    time = obj[0]
    CC = parameters[2];
    root = time.map(t => Math.exp(-t/CC));

    return root
}
//out = coneCoeffPDEVec(parameters, obj)
// console.log(coneCoeffPDEVec(parameters, obj))


function coneTauPDE(parameters, obj){
    time = obj[0]
    CC= parameters[1];
    tau = parameters[2];

    root = time.map(t => CC* Math.exp(-t/tau) *(t/(tau *tau)));
    resid = residuals(parameters,obj)

    out = -2*dot(root,resid)
    return out
}

function coneTauPDEVec(parameters, obj){
    time = obj[0]
    CC = parameters[1]
    tau = parameters[2]
    root = time.map(t => CC* Math.exp(-t/tau) *(t/(tau *tau)));

    return root
}
// out = coneTauPDEVec(parameters, obj)
// console.log(coneTauPDEVec(parameters, obj))

    

function rodS2PDE(parameters, obj){
    time = obj[0]
    alpha = parameters[4];
    b = residuals(parameters,obj);

    root = time.map(t => ((t - alpha)/(1 + Math.exp(-20*(t - alpha)))));
    out = -2*dot(root,b)
    return out
}


function rodS2PDEVec(parameters, obj){
    time = obj[0]
    alpha = parameters[4];

    root = time.map(t => ((t - alpha)/(1 + Math.exp(-20*(t - alpha)))));
   
    return root
}
//out =rodS2PDEVec(parameters, obj)
//console.log(rodS2PDEVec(parameters, obj))

function rodAlphaPDE(parameters, obj){
    time = obj[0];
    b = residuals(parameters,obj);

    S2 = parameters[3]
    alpha = parameters[4]

    root = time.map(t => (S2/(1 + Math.exp(-20 * (t - alpha))) - S2 * (t - alpha) * 
        (Math.exp(-20 * (t - alpha)) * (-20))/(1 + Math.exp(-20 * (t - alpha)))**2)
 );
  out = 2*dot(root,b)
    return out

}

function rodAlphaPDEVec(parameters, obj){
    time = obj[0];
    S2 = parameters[3]
    alpha = parameters[4]

    root = time.map(x => -(S2 / (1 + Math.exp(-20 * (x - alpha))) - S2 * (x - alpha) * 
        (Math.exp(-20 * (x - alpha)) * (-20)) / (1 + Math.exp(-20 * (x - alpha))) ** 2));
    
    return root

}

// console.log(rodAlphaPDE(parameters, obj))
//out = rodAlphaPDEVec(parameters, obj)
//console.log(rodAlphaPDEVec(parameters, obj))
function pointGrad(parameters, obj){
    x = obj[0][obj[0].length -1];
    y = obj[1][obj[1].length -1];
    p0 = parameters
    p1 = parameters.map(p => p*0.98)
    y0 = y - (p0[0] +p0[1]*Math.exp(-x/p0[2]) + p0[3]*(x - p0[4])/(1 + Math.exp(-20*(x = p0[4]))));
    y1 = y - (p1[0] +p1[1]*Math.exp(-x/p1[2]) + p1[3]*(x - p1[4])/(1 + Math.exp(-20*(x = p1[4]))));
    delY = y1*y1  - y0*y0
    delP = p0.map(t => delY/(0.02*t))
    return delP
}

//console.log(pointGrad(parameters, obj))

//build the matrix G see notes

function GradVec(parameters, obj){
    var ct= coneThresholdPDEVec(parameters, obj);
    
    var cc= coneCoeffPDEVec(parameters, obj);
   
    var tau = coneTauPDEVec(parameters, obj);
    
    var S2 = rodS2PDEVec(parameters, obj);
    
    var alpha = rodAlphaPDEVec(parameters, obj);
   G = [ct, cc, tau, S2, alpha]; 
   tGG= math.multiply( G,math.transpose(G))// js is row based, unlike R which is column based
   //OUT = math.diag(math.inv(OUT))
    return tGG
}
// console.log(GradVec(parameters, obj))



function Grad(parameters, obj){
    var ct= coneThresholdPDE(parameters, obj);
    
    var cc= coneCoeffPDE(parameters, obj);
   
    var tau = coneTauPDE (parameters, obj);
    
    var S2 = rodS2PDE(parameters, obj);
    
    var alpha = rodAlphaPDE(parameters, obj);
   
    return [ct, cc, tau, S2, alpha]; 
}
// out = math.diag(GradVec(parameters, obj))

// console.log(Grad(parameters,obj))



function errScore(parameters, obj){
    nobs = obj[0].length;
    dfP5 = nobs - 5;

    let val = leastSquare(parameters, obj);
    err = val/dfP5
    Z = critValue(dfP5);

    
    return [err, Z]
}

//console.log(errScore(parameters, obj))
//out = errScore(parameters, obj)
// console.log(obj[0])

function standardErrors(parameters, obj){
    var Variance = errScore(parameters,obj)[0]
    tGG = GradVec(parameters,obj);
    tGG = math.diag(math.inv(tGG))
    out = tGG.map(x => Variance*x)
    out = math.sqrt(out)
    return out
}
//out = standardErrors(parameters, obj)
//console.log(out)

//// SECOND DERIVATIVES AND HESSIAN

function Hessian(parameters,obj){
    
    var ct = parameters[0];
    var cc = parameters[1];
    var tau = parameters[2];
    var S2 = parameters[3];
    var alpha = parameters[4];
    
    
    var time = obj[0];
    var threshold = obj[1]
    var l11 = 2
  
    var l12 = time.map(t => 2 * Math.exp(-t/tau)).reduce((m,n)=>m+n, 0);
    var l13 = time.map(t => 2 * (cc * (Math.exp(-t/tau) * (t/tau**2)))).reduce((m,n)=>m+n, 0);
    var l14 = time.map(t => 2 * ((t - alpha)/(1 + Math.exp(-20 * (t - alpha))))).reduce((m,n)=>m+n, 0);
    var l15 = time.map(t => -(2 * (S2/(1 + Math.exp(-20 * (t - alpha))) + S2 * (t - alpha) * (Math.exp(-20 * (t - alpha)) * 20)/(1 + Math.exp(-20 * (t - alpha)))**2))).reduce((m,n)=>m+n, 0);

    var l22 = time.map(t => 2 * (Math.exp(-t/tau) * Math.exp(-t/tau))).reduce((m,n)=>m+n,0);
    var l23 = time.map((e,i) => -(2 * (Math.exp(-e/tau) * (e/(tau*tau)) * (threshold[i] - (ct + cc * Math.exp(-e/tau) + S2 * (e - alpha)/(1 + Math.exp(-20 * (e - alpha))))) - Math.exp(-e/tau) * (cc * (Math.exp(-e/tau) * (e/(tau*tau))))))).reduce((m,n)=>m+n, 0);
    var l24 = time.map(t => 2 * (Math.exp(-t/tau) * ((t - alpha)/(1 + Math.exp(-20 * (t - alpha))))) ).reduce((m,n)=>m+n, 0);
    var l25 = time.map(t => -(2 * (Math.exp(-t/tau) * (S2/(1 + Math.exp(-20 * (t - alpha))) + S2 * (t - alpha) * (Math.exp(-20 * (t - alpha)) * 20)/(1 + Math.exp(-20 * (t - alpha)))**2)))).reduce((m,n)=> m+n, 0);

    var l33 = time.map((t,i) => -(2 * (cc * (Math.exp(-t/tau) * (t/(tau*tau)) * (t/(tau*tau)) - Math.exp(-t/tau) * (t * (2 * tau)/((tau*tau))**2)) * (threshold[i] - (ct + cc * Math.exp(-t/tau) + S2 * (t - alpha)/(1 + Math.exp(-20 * (t - alpha))))) - cc * (Math.exp(-t/tau) * (t/(tau*tau))) * (cc * (Math.exp(-t/tau) * (t/(tau*tau))))))).reduce((m,n)=> m+n, 0);
    var l34 = time.map(t => 2 * (cc * (Math.exp(-t/tau) * (t/tau*tau)) * ((t - alpha)/(1 + Math.exp(-20 * (t - alpha)))))).reduce((m,n) => m+n, 0)
    var l35 = time.map(t => -(2 * (cc * (Math.exp(-t/tau) * (t/(tau*tau))) * (S2/(1 + Math.exp(-20 * (t - alpha))) + S2 * (t - alpha) * (Math.exp(-20 * (t - alpha)) * 20)/(1 + Math.exp(-20 * (t - alpha)))**2)))).reduce((m,n) => m+n, 0)

    var l44 = time.map(t=> 2 * ((t - alpha)/(1 + Math.exp(-20 * (t - alpha))) * ((t - alpha)/(1 + Math.exp(-20 * (t - alpha))))) ).reduce((m,n)=> m+n, 0);
    var l45 = time.map((t,i)=> -(2 * ((t - alpha)/(1 + Math.exp(-20 * (t - alpha))) * (S2/(1 +  Math.exp(-20 * (t - alpha))) + S2 * (t - alpha) * (Math.exp(-20 * (t - alpha)) * 20)/(1 + Math.exp(-20 * (t - alpha)))^2) - (1/(1 + Math.exp(-20 * (t - alpha))) + (t - alpha) * (Math.exp(-20 * (t - alpha)) * 20)/(1 + Math.exp(-20 * (t - alpha)))^2) * (threshold[i] - (ct + cc * Math.exp(-t/tau) + S2 * (t - alpha)/(1 + Math.exp(-20 * (t - alpha)))))))).reduce((m,n)=> m+n, 0);

    var l55= time.map((t,i)=> 2 * (((S2 * (t - alpha) * (Math.exp(-20 * (t - alpha)) * 20 * 20) - S2 * (Math.exp(-20 * (t - alpha)) * 20))/(1 + Math.exp(-20 * (t - alpha)))**2 - S2 * (t - alpha) * (Math.exp(-20 * (t - alpha)) * 20) * (2 * (Math.exp(-20 * (t - alpha)) * 20 * (1 + Math.exp(-20 * (t - alpha)))))/((1 + Math.exp(-20 * (t - alpha)))**2)**2 - S2 * (Math.exp(-20 * (t - alpha)) * 20)/(1 + Math.exp(-20 * (t -  alpha)))**2) * (threshold[i] - (ct + cc * Math.exp(-t/tau) + S2 * (t - alpha)/(1 + Math.exp(-20 * (t - alpha))))) + (S2/(1 + Math.exp(-20 * (t - alpha))) + S2 * (t - alpha) * (Math.exp(-20 * (t - alpha)) * 20)/(1 + Math.exp(-20 * (t - alpha)))**2) * (S2/(1 + Math.exp(-20 * (t - alpha))) + S2 * (t - alpha) * (Math.exp(-20 * (t - alpha)) * 20)/(1 + Math.exp(-20 * (t - alpha)))**2))
    ).reduce((m,n)=> m+n, 0);

    var r1 = [l11, l12, l13,l14,l15 ];
    var r2 = [l12, l22, l23,l24,l25 ];
    var r3 = [l13, l23, l33,l34,l35 ];
    var r4 = [l14, l24, l34,l44,l45 ];
    var r5 = [l15, l25, l35,l45,l55 ];

        
    var Hess = [r1, r2, r3,r4,r5]
    return Hess

}

//console.log(Hessian(parameters, obj))

function varianceAt_X(parameters,obj){
    Curl = pointGrad(parameters, obj)
    Var = math.inv(Hessian(parameters, obj));
    varBoost = errScore(parameters,obj)
    ptVar = math.multiply(math.multiply(Curl,Var),math.transpose(Curl))*varBoost[0]

return ptVar

}
//console.log(varianceAt_X(parameters,obj))

function intervalAt_x(parameters, obj){
    Curl = pointGrad(parameters, obj)
    Var = math.inv(Hessian(parameters, obj));
    varBoost = errScore(parameters,obj)
    Zscore = varBoost[1]
    ptVar = math.multiply(math.transpose(Curl),math.multiply(Curl,Var))*varBoost[0]
    
    Y = fiveParameterModel(parameters,obj)
    Y = Y[Y.length - 1]
    ub = Y + Zscore * Math.sqrt(ptVar)
    lb = Y - Zscore * Math.sqrt(ptVar)

    return [lb, ub]
}

//console.log(intervalAt_x(parameters, obj))


function parameterCIs(parameters, obj){
    tNumber = errScore(parameters,obj)[1]
    standErr = standardErrors(parameters, obj)
    band = standErr.map(s => tNumber*s);
    
    ub = parameters.map((p,i) => p+band[i]);
    lb = parameters.map((p,i) => p-band[i]);
    truth = ub.map((a,i) => a/lb[i]>0?1:0)
    ub = truth.map((e,i) => e == 0?+2:ub[i])
    lb = truth.map((e,i) => e == 0?-2:lb[i])
    return [lb,ub]
}
//ULUB = parameterCIs(parameters, obj);
//console.log(ULUB)
//out = ULUB

fn = "Math.pow((threshold[i] - (ct + cc * exp(-time/tau) + S2 * (time - alpha)/(1 + exp(-20* (time - alpha))))),2)"

fnc = "Math.pow((threshold - (ct + cc * exp(-time/tau))),2)"

function setCITable(parameters, obj){
return null
}
