(* ::Package:: *)

ClearAll["Global`*"]

X = {x1, x2, x3};
Y = {y1, y2, y3};

norm[v_]:=Sqrt[v . v]

A0 = -k3^2(kt^2-k1^2)/((kt^2-kc^2)kt^2);
A1 = -k3^2(kc^2-k1^2)/((kc^2-kt^2)kc^2);

g[x_,y_,k_] := Exp[I k norm[x-y]] / (4Pi norm[x-y])
g\[Theta][R_,k_] := Exp[I k R] / (4Pi R)

u[i_,j_, x_, y_, \[Mu]_] := 1/(\[Mu] ks^2) (D[A0 g[x,y,kt] + A1 g[x,y, kc] + g[x,y,ks] , x[[i]], x[[j]]] + ks^2 KroneckerDelta[i, j] g[x,y, ks])
\[Theta][i_, x_,y_,\[Mu]_, \[Lambda]_] := B / A (kt^2-kc^2) D[g[x,y,kc] - g[x,y,kt], x[[i]]];
\[Theta]n[i_, x_, y_,\[Mu]_, \[Lambda]_] := B / A (kt^2-kc^2)  D[D[g\[Theta][R,kc] - g\[Theta][R,kt], R] /. R-> norm[x-y], x[[i]]];

\[Sigma][i_,j_, k_,x_,y_, \[Mu]_, \[Lambda]_]  := \[Lambda] KroneckerDelta[j, k] (D[u[i, 1, x, y, \[Mu]], x[[1]]] +D[u[i, 2, x, y, \[Mu]], x[[2]]] + D[u[i, 3, x, y, \[Mu]], x[[3]]]) + \[Mu] (D[u[i, j, x, y, \[Mu]], x[[k]]] +D[u[i, k, x, y, \[Mu]], x[[j]]]) - B \[Theta][i,x,y,\[Mu], \[Lambda]] KroneckerDelta[j, k]

node[n_, m_, r_] := {r Sin[Pi n/nmax]Cos[2 Pi m/mmax], r Sin[Pi n/nmax]Sin[2 Pi m/mmax], r Cos[Pi n/nmax]}
quad[n_,m_,r_] := {node[n,m,r], node[n, m+1,r], node[n+1, m,r], node[n+1,m+1,r]}

y[\[Xi]_, \[Eta]_, yi_] := 1/4(1-\[Xi])(1-\[Eta]) yi[[1]] + 1/4(1+\[Xi])(1-\[Eta]) yi[[2]] + 1/4(1-\[Xi])(1+\[Eta]) yi[[3]] +1/4(1+\[Xi])(1+\[Eta]) yi[[4]]

(* rigid *)
f[i_, j_] := u[i, j, X, Y, mu] *norm[Cross[D[y[xi, eta, quad[nq,mq,1]],xi],  D[y[xi, eta, quad[nq,mq,1]],eta]]]
f2[i_] := \[Theta][i, X, Y, mu, lambda]*norm[Cross[D[y[xi, eta, quad[nq,mq,1]],xi],  D[y[xi, eta, quad[nq,mq,1]],eta]]]

(* temp BIR *)
\[Phi]s[R_] := k1^2/ cv B /(A (kt^2-kc^2))*(g\[Theta][R, kt]-g\[Theta][R, kc])
\[Theta]s[R_] := 1/((kt^2-kc^2))*((kt^2-k3^2)g\[Theta][R, kt]-(kc^2-k3^2) g\[Theta][R, kc])

temp := (\[Theta]s[R] /. R -> Norm[X-Y]// FullSimplify)norm[Cross[D[y[xi, eta, quad[nq,mq,1]],xi],  D[y[xi, eta, quad[nq,mq,1]],eta]]]
disp[i_] := (Grad[Evaluate[\[Phi]s[R] /. R-> norm[X-Y]], {x1,x2,x3}]norm[Cross[D[y[xi, eta, quad[nq,mq,1]],xi],  D[y[xi, eta, quad[nq,mq,1]],eta]]] // Simplify)[[i]]

(* Using the ToMatlab package by Harri Ojanen (1999) https://library.wolfram.com/infocenter/MathSource/577/
the functions f[i,j] for i,j={1,2,3} form the entries of Uij_vector_BIR.m,
f2[i] for i = {1,2,3} are the entries of temp_vector_BIR.m, and
disp[i] for i = {1,2,3} and temp are the entries of scalar_rigid_BIR.m  *)
