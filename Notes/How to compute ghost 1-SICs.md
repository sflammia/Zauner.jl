# How to compute ghost 1-SICs.

### What defines a ghost?

Let $d = d_h >3$ be an integer. Then define $\Delta_0$ to be the fundamental discriminant of $\mathbb{Q}\bigl(\sqrt{(d+1)(d-3)}\bigr).$ From the factorization $(d+1)(d-3) = f_h^2 \Delta_0$, we implicitly define the positive integer $f_h$. Then a ghost can be specified by giving

- An integer $d = d_h > 3$. 
- A quadratic form $Q$ with discriminant $\Delta_Q = f_Q^2\Delta_0$, where the conductor $f_Q$ divides $f_h$.

Remarks: The ghost SIC formula depends on $Q$ (via the stabilizer of $Q$), however any $Q$ in the same $\mathrm{SL}_2(\mathbb{Z})$ class give ghosts on the same extended Clifford orbit. Therefore, we can choose $Q$ to be reduced wlog. This is motivated by the fact that non-reduced forms require greater computational resources in general.

To define a SIC (rather than just a ghost), we also need to specify a sign-switching Galois automorphism $g(\Delta_0) = -\Delta_0$. This could be specified by an explicit choice of minimal polynomial for the field containing the SIC as well as an interval containing a root. 

One additional strategy for minimizing the computational effort is to work not just with a reduced form $Q$, but with one whose stabilizer has a decomposition into $S$ and $T$ generators that is as short as possible, with nonnegative powers of $T$ only (except for maybe on the end). Thus, one should choose such a _minimized_ reduced form $Q$. Even this is not unique, not even up to cycling around a pure period, so a choice has to be made. 

### How do I compute a ghost?

There are many steps required in order to compute the ghost. The first step is to compute $L = L(Q)$, a generator of the stability group of $Q$, i.e. the group of matrices $M \in \mathrm{SL}_2(\mathbb{Z})$ such that $M^T Q M = Q$. If we can write an arbitrary $M$ in this group as $M =\pm L^k$ for $k\in \mathbb{Z}$, then $L$ is a valid generator and any choice (there are four) gives us a choice of $L(Q)$. (The freedom in choice is $\pm L,\pm L^{-1}$.)

To compute $L(Q)$, we choose a fundamental unit of the order with discriminant $\Delta_Q$ and conductor $f_Q$. In PARI/GP, this is done with `uQ=quadunit(DeltaQ)` to get a completely positive unit with $u_Q>1$. If the norm is negative, then square it to get a unit with positive norm, $u_Q \leftarrow  u_Q^2$. Call this (possibly squared) unit $u_Q$. This unit can be split into its rational and irrational components, but a more convenient basis is the "canonical" integral basis $(1,\omega)$, where $\omega = \frac{(\Delta_Q \bmod 4)+\sqrt{\Delta_Q}}{2}$.  In PARI/GP, the components in this basis are given by `real` and `imag`. Setting $x = \mathrm{tr}(u_Q)/2$ and `y=imag(uQ)`, then an explicit formula for a stabilizer generator is $L(Q) = x I + y SQ$, where $S = [0,-1;1,0]$ is a generator of  $\mathrm{SL}_2(\mathbb{Z})$ and $Q$ here is understood to be the matrix of the quadratic form $Q$. 

The height of the tower, $h$, is given as follows. Let $u_0$ be the fundamental unit for $\mathbb{Z}(\Delta_0)$, and as before if the norm is negative, then square it to get a unit with positive norm, $u_0 \leftarrow  u_0^2$. Because of the equation $d_h = 1+u_0^h + u_0^{-h}$, we have $h = \mathrm{acosh}((d_h - 1)/2)/\log(u_0)$. 

Next, compute $A = L^{3h/h_0}$. This is the smallest power of $L$ such that $A \equiv I \bmod d$. Here $h_0$ is the smallest height such that $f_Q | f_{h_0}$. That is, for all $k$, we have $(d_k+1)(d_k-3) = f_k^2 \Delta_0$. For some least value of $k$ in that equation, call it $h_0$, we have that $f_Q |  f_{h_0}$. It might be computationally faster to find the order of $L \bmod d$. 

Decompose $A = T^{e_1}ST^{e_2}\cdots ST^{e_{n+1}}$. If we made a sensible choice of sign in our choice of $Q$ and $L(Q)$, then there is no possibility of $-I$ appearing in this decomposition. This decomposition needs to have $e_k \ge 2$ for all $1 < k \le n$ as a sufficient condition for avoiding singular behavior of subsequent functions. This can be done by repeated application of the GCD algorithm using the largest nonnegative integer (instead of the nearest integer). 

Next, compute $\psi(A)$, the Rademacher invariant of $A$, which is a class invariant of $\mathrm{SL}_2(\mathbb{Z})$. (The closely related Meyer invariant $\phi(A) = -\psi(A)/3$ for hyperbolic elements $A$, and is given by Eq. (5.34) in draft.tex.) An explicit form for an element of $\gamma \in \mathrm{SL}_2(\mathbb{Z})$ with $\gamma = \begin{pmatrix}a&b\\c&d\end{pmatrix}$ is given as follows. Let $t = \mathrm{tr(\gamma)}$ be the matrix trace of $\gamma$, so $t = a+d$. Then 

###### (**COMMENT**:  I think this may assume $c>0$ and $t>2$.  STF:Probably not c > 0 since there is an explicit sign(c) in one of the equations. But t > 2 maybe. I got this from Ghys (I think). )

​    $\displaystyle \psi(\gamma) = \begin{cases}b/d & \mbox{if }c=0,\\ t/c - 12\,
\mathrm{sign}(c) s(a,|c|) - 3\,\mathrm{sign}(ct) & \mbox{otherwise}.\end{cases} $ 

Here $s(a,b)$ is the Dedekind sum, 

​    $\displaystyle s(a,b) = \sum_{n=1}^{b-1} \Bigl(\!\Bigl(\frac{n}{b}\Bigr)\!\Bigr) \Bigl(\!\Bigl(\frac{an}{b}\Bigr)\!\Bigr),$ 

where $(\!(x)\!) = 0$ if $x \in \mathbb{Z}$ and $(\!(x)\!) = x-\lfloor x\rfloor-1/2$ otherwise. 

Define $\beta$ to be a root of $[1,\beta] Q [1,\beta]^T = 0$.  Either choice works, but we choose the root > 1. That is, if $Q = \langle a,b,c\rangle$, then $\beta = (-b + \sqrt{b^2-4ac})/2a$ or equivalently

​    $$\displaystyle \beta = \frac{-b + \sqrt{\Delta_Q}}{2a}\,.$$ 

Now define a symplectic inner product as follows: $\langle\!\langle \boldsymbol{p},\boldsymbol{q} \rangle\!\rangle := p_2 q_1 - p_1 q_2$. Furthermore, we define $\langle\!\langle\boldsymbol{p},\beta\rangle\!\rangle = p_2 \beta - p_1$. 

We next define the "exponential" q-Pochhammer symbol as follows. First define $e(z) = \exp(2\pi i z)$. Then we have: 

​    $$\displaystyle \varpi(z,\beta) = \prod_{j=0}^\infty (1-e(z+j\beta))$$ 

and for integers $m$ we also define:

​    $$\displaystyle \varpi_m(z,\beta) = \begin{cases}
\prod_{j=0}^{m-1}(1-e(z+j\beta)) & m \ge 1\\ 
1 & m=0\\ 
\prod_{j=m}^{-1} (1-e(z+j\beta))^{-1} & m \le -1
\end{cases}\,.$$

Now define $j_L(\beta) = c\beta+d$ for $L = \begin{pmatrix}a&b\\ c&d\end{pmatrix}$. We also define $\displaystyle L\beta = \frac{a\beta + b}{c\beta+d}$. 

Next, define

​    $$\displaystyle \sigma_L(z,\beta) = \frac{\varpi\bigl(\frac{z}{j_L(\beta)},L\beta\bigr)}{\varpi(z,\beta)}.$$

Now define $\rho$, aka "shin", but I don't think I can typeset shin in markdown. :)  

​    $$\displaystyle \rho_{A,d}(\boldsymbol{p},\beta) = \frac{\sigma_A(\langle\!\langle A^{-1}\boldsymbol{p},\beta \rangle\!\rangle/d,\beta )}{\varpi_m(\langle\!\langle \boldsymbol{p},\beta \rangle\!\rangle/d,\beta)}$$ 

Where $m$ is given in terms of the elements of $A$ as $m = (-c p_1+(a-1)p_2)/d$. Next define:

​    $$\displaystyle s_{\boldsymbol{q}}(\boldsymbol{p}) = 
\begin{cases}
1 & d \mbox{ odd}\\ 
(1+p_1)(1+p_2) + \langle\boldsymbol{q},\boldsymbol{p}\rangle & d \mbox{ even}
\end{cases}$$ 

and if $A = \begin{pmatrix}a&b\\ c&d\end{pmatrix}$, define

​    $$\displaystyle Q_A = \frac{1}{d_h (d_h-2)}\langle c, d-a, -b\rangle = \frac{1}{d_h (d_h-2)}\begin{pmatrix}
c & (d-a)/2\\
(d-a)/2 & -b
\end{pmatrix}$$,

and define $\xi_d = -\exp(\pi i/d)$. For $\boldsymbol{q}$ one of $[0,0], [0,1], [1,0], [1,1]$, we finally define the following root of unity

​    $$\displaystyle \omega(A,d_h,\boldsymbol{p},\boldsymbol{q}) = (-1)^{s_{\boldsymbol{q}}(\boldsymbol{p})} e(-\psi(A)/24) \xi_d^{-\boldsymbol{p}^T Q_A \boldsymbol{p}}\,.$$

The ghost overlap function is then given by

​    $$\displaystyle \nu_{A,\beta,\boldsymbol{q},d}(\boldsymbol{p}) = \begin{cases}
1 & \boldsymbol{p} = \boldsymbol{0} \bmod d\\
\frac{\omega(A,d,\boldsymbol{p},\boldsymbol{q})}{\sqrt{d+1}} \rho_{A,d}(\boldsymbol{p},\beta) & \mbox{otherwise}\,.
\end{cases}$$ 

With these definitions, we can construct the ghost projector as

​    $$\displaystyle \Pi_{A,\beta,\boldsymbol{q},d} = \frac{1}{d} \sum_{\boldsymbol{p}} \nu_{A,\beta,\boldsymbol{q},d}(\boldsymbol{p}) D_{\boldsymbol{p}}$$

Where $D_{\boldsymbol{p}}$ is the displacement operator. 

To expedite the computation of $\sigma_A$, we can use the following product theorem. Consider that $A = T^{e_1}ST^{e_2}\cdots ST^{e_{n+1}}$ as before, and we can always choose $e_k \ge 2$ for $1 < k < n+1$. For our choices of $A$ above (inherited from our choice of $L(Q)$), we never have to worry about the case $A = -T^e$.  Define $A_k = T^{e_k}S\ldots ST^{e_{n+1}}$ to be the partial products that form $A$. Then the factorization theorem states that 

​    $$\displaystyle \sigma_A(z,\beta) = \prod_{k=2}^{n+1} \sigma_S\biggl(\frac{z}{j_{A_k}(\beta)},A_k\beta\biggr)\,,$$

where $S = \begin{pmatrix}0&-1\\ 1&0\end{pmatrix}$ is a standard generator for $\mathrm{SL}(2,\mathbb{Z})$. A nice way to write this is to define $z_k = \frac{z}{j_{A_k}(\beta)}$ and $\beta_k = A_k\beta$, so that the product looks like

​    $$\displaystyle \sigma_A(z,\beta) = \prod_{k=2}^{n+1} \sigma_S(z_k,\beta_k)\,,$$

Then the problem reduces to computing a set of functions of the form $\sigma_S(z_k,\beta_k)$. Subsequently, this can be done using the formula

​    $$ \displaystyle \sigma_S(z,\beta) = e\biggl(\frac{6z^2+6(1-\beta)z+\beta^2-3\beta+1}{24\beta}\biggr) S_2(z+1,\beta)\,.$$

Here $S_2(z,(\omega_1,\omega_2))$ is the double sine function, and $S_2(z,\beta) = S_2(z,(\beta,1))$ is the "reduced" double sine function where the second period is set to 1. 

When $S_2(z,(\omega_1,\omega_2))$ takes arguments obeying $0<\mathrm{Re}(z) < \mathrm{Re}(\omega_1+\omega_2)$, then we have the following integral representation:

​    $$\displaystyle S_2(z,(\omega_1,\omega_2)) = \exp\biggl(\int_0^\infty \biggl(\frac{\sinh\bigl((\omega_1+\omega_2-2z)t/2\bigr)}{2\sinh(\omega_1 t/2)\sinh(\omega_2 t/2)} - \frac{\omega_1+\omega_2 - 2 z}{t} \biggr) \frac{\mathrm{d}t}{t} \biggr) \,.$$  

Note that this definition comes from Shintani, and Kurokawa and others would use a minus sign in front of the integral. We should adjust our formulas accordingly to match the contemporary (Kurokawa) definitions.  When the double sine function does not satisfy that domain constraint, then we can use the transformation laws to shift the double sine into that domain, picking up additional factors. Specifically, we have the following transformation law.

​    $$ \displaystyle \sigma_S(z,\beta) = \sigma_S(z+m\beta+n,\beta)\varpi_{-n}\biggl(\frac{z}{\beta},\frac{-1}{\beta}\biggr) \varpi_{m}^{-1}(z,\beta)$$

Since $\beta > 0$, the correct domain for us to shift to can always be found with the following choices: $m = 0, n=\lfloor -z\rfloor$. In that case, we have 

​    $$ \displaystyle \sigma_S(z,\beta) = \sigma_S\bigl(z+\lfloor -z \rfloor,\beta\bigr)\,
\varpi_{-\lfloor -z \rfloor}\biggl(\frac{z}{\beta},\frac{-1}{\beta}\biggr)$$



### Hirzebruch-Jung reduced forms

For a quadratic form $Q = ⟨ a,b,c\rangle$ with discriminant $\Delta$ and $L \in \mathrm{GL}(2,\mathbb{Z})$, we define
    $$ Q_L = (\det L) L^T Q L \ , \quad \beta_{Q,\pm} = \frac{-b\pm \sqrt{\Delta}}{2a}\,.$$
Recall of course that $\Delta$ is a class invariant, so although $\beta_{Q,\pm}$ changes under the action of $L$, only $b$ and $a$ change. 

A result that is useful for us is that $\beta_{Q,\pm}$ has an HJ-continued fraction expansion that is purely periodic if and only if $\beta > 1 > \beta^\prime > 0$, where $\beta^\prime$ is the Galois conjugate of $\beta$. 

Recall the notion of ``ordinary'' or Euclidean reduction. A quadratic form $Q$ is Euclidean reduced if 
    $$0 < \sqrt{\Delta} − b < 2|a| < \sqrt{\Delta} + b\,.$$
The number $b+\sqrt{\Delta}/2|a|$ has a purely periodic Euclidean continued fraction expansion if and only if $Q$ is Euclidean reduced. Corresponding to this we say $Q$ is Hirzebruch-Jung (HJ-) reduced if
    $$0 < −\sqrt{\Delta} − b < 2|a| < \sqrt{\Delta} − b\,.$$
The number $−b+\sqrt{\Delta}/2|a|$ has a purely periodic HJ-continued fraction expansion if and only if $Q$ is
HJ-reduced. 

Now notice that $Q_J := \langle -a,b,-c\rangle$ is conjugate to $Q$. Morefover, it is HJ reduced iff $Q$ is HJ reduced. Therefore we can consider $a > 0$ without loss of generality. Then the map taking $Q$ to $\beta_{Q,+}$ is a bijection between the set of HJ reduced forms and purely periodic HJ continued fractions. 

Now let $W$ be the set of forms on a given $\mathrm{GL}(2,\mathbb{Z})$ orbit such that $a > 0$. Then the following two sets are precisely the sets of Euclidean or HJ reduced forms in $W$. 
    $$ W_E = \{Q\in W : \beta_{Q,-} < -1 < 0 < \beta_{Q,+} < 1\}$$
    $$ W_{HJ} = \{Q\in W : 0< \beta_{Q,-} < 1 < \beta_{Q,+} < 1\}$$
Then $W_E$ and $W_{HJ}$ are precisely the set of Euclidean and HJ reduced forms respectively in $W$. 





###### (Insert transformation laws here to make sure that we agree with Appleby.)

### Comments:

This exercise has illuminated for me some weaknesses that we want to address when we write up the paper. First, we don't have anything close to even this meager version of a straight-line view of the entire construction. This is essential. Second, It seems to me that many of these moving parts from the construction could be profitably synthesized into tighter units if only we could manage to get the entire thing into view to see what best to combine. So I hope that progress on the first point will lead to a simpler form of the construction. Third, I've said this before, but we really should somehow make explicit that the ghost transforms covariantly with the Clifford group, so that the ghosts modulo extended Clifford group elements are class invariants. We've proven one direction: that if two points are SL(2,Z) related, then they are also Clifford related. To prove the converse, maybe we can use the triple products? Fourth, I know enough by now to know that there is quite a lot of geometry happening here under the surface. Perhaps this is too much to ask for the present paper, but I feel like this geometric perspective needs a much more prominent role in our discussions. Personally, I believe that understanding the geometric perspective is the only way that I'll ever really intuitively grasp this construction. 
