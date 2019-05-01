##<center>Assignment 3</center>
<center>成子浩 2018212327</center>

<font size=4>**Problem 1 Quantile Regression**</font>

a) *QP*
Firstly, we could introduce a some slack variables \(\xi_i \) to replace the absolute value loss. Then the problem can be formulated as 
$$
min\quad \frac{\lambda}{2}||w||^2 + \frac{1}{n}\sum_{i=1}^n\xi_i\\ s.t.\quad |y_i - w^Tx_i - b| = \xi_i\\ \xi_i \geq 0 
$$
&ensp;&ensp;That is not a QP since the constraint is not linear. However, we found that the problem is equivalent to the following formulation. 
$$
min\quad \frac{\lambda}{2}||w||^2 + \frac{1}{n}\sum_{i=1}^n\xi_i\\ 
s.t.
\begin{aligned}
\quad y_i - w^Tx_i - b &\leq \xi_i\\
\quad y_i - w^Tx_i - b &\geq -\xi_i\\
\quad \xi_i &\geq 0
\end{aligned}
$$

&ensp;&ensp;For a certain \(i_{th}\) sample point, it has to satisfy either the first inequality or the second one, or the objective value doesn't reach its minimum state. Thus the quadratic programing is derived.

b)*Dual form*
Firstly, we introduce \(\mu\) and \(\theta\) as the dual variables. So the *Lagrangian function* can be formulated as 
$$
\begin{aligned}
min\quad L(w,b,\mu,\theta) 
&=min\quad \frac{\lambda}{2}||w||^2+\frac{1}{n}\sum_{i=1}^n\xi_i+\sum_{i=1}^n\mu_i(y_i-w^Tx_i-b-\xi_i)-\sum_{i=1}^n\theta_i(y_i-w^Tx_i-b+\xi_i)\\
&=min\quad \sum_{i=1}^n(\mu_i-\theta_i)(y_i-w^Tx_i-b)+\frac{\lambda}{2}||w||^2+\sum_{i=1}^n(\frac{1}{n}-\mu_i-\theta_i)\xi_i\\
&=\begin{cases}\sum_{i=1}^n(\mu_i-\theta_i)(y_i-w^Tx_i-b)+\frac{\lambda}{2}||w||^2\quad \quad if\quad \frac{1}{n}\geq\mu_i+\theta_i\\
-\infty \quad\quad o.w.
\end{cases}
\end{aligned}
$$
&ensp;&ensp;By soluting \(\frac{\partial L(w,b,\mu,\theta)}{\partial w} = 0\) and \(\frac{\partial L(w,b,\mu,\theta)}{\partial b} = 0\) , we can derive following equations, 
$$
\begin{aligned}
w &= \frac{1}{\lambda}\sum_{i=1}^n(\mu_i-\theta_i)x_i\\
0 &= \sum_{i=1}^n(\mu_i-\theta_i)
\end{aligned}
$$
&ensp;&ensp;Then we plug \(w^*\) back into the dual objective to element \(w\) and \(b\), i.e. 
$$
min\quad L(w,b,\mu,\theta) = \sum_{i=1}^n(\mu_i-\theta_i)y_i-(1-\frac{1}{2\lambda})\sum_{i=1}^n\sum_{j=1}^n(\mu_i-\theta_i)(\mu_j-\theta_j)x_i^Tx_j
$$
&ensp;&ensp;So the dual form of the primal problem is given by
$$
Dual:max\quad \sum_{i=1}^n(\mu_i-\theta_i)y_i-(1-\frac{1}{2\lambda})\sum_{i=1}^n\sum_{j=1}^n(\mu_i-\theta_i)(\mu_j-\theta_j)x_i^Tx_j\\s.t.\quad\quad\quad
\begin{aligned}
\mu_i+\theta_i &\leq \frac{1}{n} \quad \forall i\\
\sum_{i=1}^n(\mu_i-\theta_i)&=0\\
\mu_i,\theta_i &\geq 0\quad \forall i
\end{aligned}
$$

c)*Primal solution*
From *b)* we have derived the formulation of \(w^*\), i.e. \(w^* = \frac{1}{\lambda}\sum_{i=1}^n(\mu_i-\theta_i)x_i\). Firstly, according to **KKT conditions**, if there exists a \(\mu_i\) satisfying \(\mu_i \not= 0\), then we have \(\xi_i = y_i - w^Tx_i - b\). Samely, if \(\theta_i \not= 0\), then \(\xi_i = -(y_i - w^Tx_i - b)\). We have claimed in *b)* that for a certain sample \((x_i,y_i)\), it either satisfies \(\mu_i \not= 0\) or \(\theta_i \not=0\). So we can use \(x_i,y_i,w,b\) to represent \(\xi_i\) for all i without **absolute value**. We define two new set as \(S_1 = \{i:\mu_i\not=0\}\) and \(S_2 = \{j:\theta_j\not=0\}\). Then we plug \(w^*\) and \(\xi_i\) back into the primal objective to get the following formulation
$$
d^* = \frac{\lambda}{2}||w||^2+\frac{1}{|S1|}\sum_{i\in S1}(y_i-w^Tx_i-b)-\frac{1}{|S2|}\sum_{j\in S2}(y_j-w^Tx_j-b),
$$
where \(d^*\) denotes the optimal value of the dual problem. Since \(p^* = d^*\), and all the other coefficients except \(b\) could be computed by dual problem or be represented by some combination of the known coefficients. We could use the above equation to compute \(b^*\).
<br>

<font size=4>**Problem 2 M-W**</font>

*a)*
The Sherman-Morrison-Woodbury formula is given by
$$
(A+uv^T)^{-1} = A^{-1} - \frac{A^{-1}uv^TA^{-1}}{1+v^TA^{-1}u},
$$
where \(A\) is an invertible square matrix and \(u,v \in \R^n\) are column vectors. In our problem, we regard \((1-\rho)Diag(\vec\sigma)^2\) as \(A\), \(\rho\sigma\) as \(u\) and \(\sigma\) as \(v\). Then we have
$$
\begin{aligned}
V^{-1}
&=\frac{1}{1-\rho}Diag(\theta)^2 - \frac{\frac{1}{1-\rho}Diag(\theta)^2\rho\sigma\sigma^T\frac{1}{1-\rho}Diag(\theta)^2}{1+\sigma^T\frac{1}{1-\rho}Diag(\theta)^2\rho\sigma}\\
&=\frac{1}{1-\rho}Diag(\theta)^2 - \frac{\rho}{(1-\rho)(1+(n-1)\rho)}\theta\theta^T
\end{aligned},
$$
where \(\theta = [\frac{1}{\sigma_1},\frac{1}{\sigma_2},\cdots,\frac{1}{\sigma_n}]^T\).

*b)*
