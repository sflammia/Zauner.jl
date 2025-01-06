# Admissible tuple data tables

This table contains algebraic data for each inequivalent admissible tuple in dimensions 4--35, comprising 100 total tuples.
This list is conjecturally complete for all Weyl-Heisenberg covariant SICs in these dimensions.
Each admissible tuple is specified by a dimension $d$ and an integer binary quadratic form $Q$ as follows.
First factorize $(d+1)(d-3)=s^2\Delta_0$ where $\Delta_0$ is a fundamental discriminant.
Then $ (d,Q) $ gives an admissible tuple if $\mathrm{disc}(Q) = f^2\Delta_0$ where $f$ divides $s$.
The other columns can be computed from these data, but they may be difficult to compute,
for example requiring integer factoring or finding a fundamental unit.
The column $\Delta_0$ contains the fundamental discriminant of $Q$ and $h$ is the order of the class group $\mathrm{Cl}(\mathcal{O}_f)$, given in the next two columns respectively.
The Galois group $\mathrm{Gal}(E_t^{(1)}/H)$ of the candidate overlap field ramified at the first infinite place over the class field $H$ of $K=\mathbb{Q}(\sqrt{\Delta_0})$ is given in the next column.
As both the class group and the Galois group are finite and abelian, we give the canonical decomposition into cyclic groups $C_k$ of order $k$.
The column $L^n$ contains a generator $L$ of the stability group of $Q$ in $\mathrm{GL}_2(\mathbb{Z})$ and its order $n$ in $\mathrm{GL}_2(\mathbb{Z}/\bar{d})$;
that is, treating $Q$ as a symmetric matrix we have $L^T Q L = \det(L) Q$ and $L^n = 1\ (\bmod\ \bar{d})$.
If the tuple has antiunitary symmetry, we denote this with a Y in the a.u. column.
Finally, $\ell$ is the length of the word expansion of $L^n$ using the Hirzebruch-Jung (negative regular) reduction into the standard
($S$ and $T$) generators of $\mathrm{SL}_2(\mathbb{Z})$.
This is one measure of the complexity of constructing the actual fiducial vector for that input.
The $Q$ in this list were selected among class representatives to minimize this complexity, although this choice is not generally unique.
The data here are sufficient to compute a ghost fiducial in each class,
but to fully specify a SIC,
one must additionally choose a sign-switching Galois automorphism $\sqrt{\Delta_0}\to-\sqrt{\Delta_0}$ over an appropriate field extension of $K$.

| $d$ | $\Delta_0$ | $f$ | $h$ | $\mathrm{Cl}(\mathcal{O}_f)$ | $\mathrm{Gal}(E/H)$ | $Q$ | $L^n$ | $\text{a.u.}$ | $\ell$ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| $4$ | $5$ | $1$ | $1$ | $ C_{1} $ | $C_{2}^{2}$  | $\langle1,-3,1\rangle$|$\left(\begin{smallmatrix}2&-1\\1&-1\end{smallmatrix}\right)^{6}$|$\text{Y}$|$4$|
| $5$ | $12$ | $1$ | $1$ | $ C_{1} $ | $C_{8}$  | $\langle1,-4,1\rangle$|$\left(\begin{smallmatrix}4&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| $6$ | $21$ | $1$ | $1$ | $ C_{1} $ | $C_{2}\times{}C_{6}$  | $\langle1,-5,1\rangle$|$\left(\begin{smallmatrix}5&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| $7$ | $8$ | $1$ | $1$ | $ C_{1} $ | $C_{6}$  | $\langle2,-4,1\rangle$|$\left(\begin{smallmatrix}3&-1\\2&-1\end{smallmatrix}\right)^{6}$|$\text{Y}$|$7$|
| | | $2$ | $1$ | $ C_{1} $ | $C_{2}\times{}C_{6}$  | $\langle1,-6,1\rangle$|$\left(\begin{smallmatrix}6&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| $8$ | $5$ | $1$ | $1$ | $ C_{1} $ | $C_{2}\times{}C_{4}$  | $\langle1,-3,1\rangle$|$\left(\begin{smallmatrix}2&-1\\1&-1\end{smallmatrix}\right)^{12}$|$\text{Y}$|$7$|
| | | $3$ | $1$ | $ C_{1} $ | $C_{2}\times{}C_{4}^{2}$  | $\langle1,-7,1\rangle$|$\left(\begin{smallmatrix}7&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| $9$ | $60$ | $1$ | $2$ | $C_{2}$  | $C_{3}\times{}C_{6}$  | $\langle1,-8,1\rangle$|$\left(\begin{smallmatrix}8&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle5,-10,2\rangle$|$\left(\begin{smallmatrix}9&-2\\5&-1\end{smallmatrix}\right)^{3}$||$7$|
| $10$ | $77$ | $1$ | $1$ | $ C_{1} $ | $C_{2}\times{}C_{24}$  | $\langle1,-9,1\rangle$|$\left(\begin{smallmatrix}9&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| $11$ | $24$ | $1$ | $1$ | $ C_{1} $ | $C_{40}$  | $\langle3,-6,1\rangle$|$\left(\begin{smallmatrix}11&-2\\6&-1\end{smallmatrix}\right)^{3}$||$7$|
| | | $2$ | $2$ | $C_{2}$  | $C_{40}$  | $\langle1,-10,1\rangle$|$\left(\begin{smallmatrix}10&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle3,-12,4\rangle$|$\left(\begin{smallmatrix}11&-4\\3&-1\end{smallmatrix}\right)^{3}$||$7$|
| $12$ | $13$ | $1$ | $1$ | $ C_{1} $ |  $C_{2}^{4}$  | $\langle3,-5,1\rangle$|$\left(\begin{smallmatrix}4&-1\\3&-1\end{smallmatrix}\right)^{6}$|$\text{Y}$|$10$|
| | | $3$ | $1$ | $ C_{1} $ | $C_{2}^{3}\times{}C_{6}$  | $\langle1,-11,1\rangle$|$\left(\begin{smallmatrix}11&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| $13$ | $140$ | $1$ | $2$ | $C_{2}$  | $C_{4}\times{}C_{12}$  | $\langle1,-12,1\rangle$|$\left(\begin{smallmatrix}12&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle7,-14,2\rangle$|$\left(\begin{smallmatrix}13&-2\\7&-1\end{smallmatrix}\right)^{3}$||$7$|
| $14$ | $165$ | $1$ | $2$ | $C_{2}$  | $C_{2}\times{}C_{6}^{2}$  | $\langle1,-13,1\rangle$|$\left(\begin{smallmatrix}13&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle5,-15,3\rangle$|$\left(\begin{smallmatrix}14&-3\\5&-1\end{smallmatrix}\right)^{3}$||$7$|
| $15$ | $12$ | $1$ | $1$ | $ C_{1} $ | $C_{24}$  | $\langle1,-4,1\rangle$|$\left(\begin{smallmatrix}4&-1\\1&0\end{smallmatrix}\right)^{6}$||$7$|
| | | $2$ | $1$ | $ C_{1} $ | $C_{2}\times{}C_{24}$  | $\langle4,-8,1\rangle$|$\left(\begin{smallmatrix}15&-2\\8&-1\end{smallmatrix}\right)^{3}$||$7$|
| | | $4$ | $2$ | $C_{2}$  | $C_{2}\times{}C_{24}$  | $\langle1,-14,1\rangle$|$\left(\begin{smallmatrix}14&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle11,-18,3\rangle$|$\left(\begin{smallmatrix}16&-3\\11&-2\end{smallmatrix}\right)^{3}$||$10$|
| $16$ | $221$ | $1$ | $2$ | $C_{2}$  | $C_{2}\times{}C_{8}^{2}$  | $\langle1,-15,1\rangle$|$\left(\begin{smallmatrix}15&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle7,-19,5\rangle$|$\left(\begin{smallmatrix}17&-5\\7&-2\end{smallmatrix}\right)^{3}$||$10$|
| $17$ | $28$ | $1$ | $1$ | $ C_{1} $ | $C_{96}$  | $\langle2,-6,1\rangle$|$\left(\begin{smallmatrix}17&-3\\6&-1\end{smallmatrix}\right)^{3}$||$7$|
| | | $3$ | $2$ | $C_{2}$  | $C_{96}$  | $\langle1,-16,1\rangle$|$\left(\begin{smallmatrix}16&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle9,-18,2\rangle$|$\left(\begin{smallmatrix}17&-2\\9&-1\end{smallmatrix}\right)^{3}$||$7$|
| $18$ | $285$ | $1$ | $2$ | $C_{2}$  | $C_{3}\times{}C_{6}^{2}$  | $\langle1,-17,1\rangle$|$\left(\begin{smallmatrix}17&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle13,-21,3\rangle$|$\left(\begin{smallmatrix}19&-3\\13&-2\end{smallmatrix}\right)^{3}$||$10$|
| $19$ | $5$ | $1$ | $1$ | $ C_{1} $ | $C_{18}$  | $\langle1,-3,1\rangle$|$\left(\begin{smallmatrix}2&-1\\1&-1\end{smallmatrix}\right)^{18}$|$\text{Y}$|$10$|
| | | $2$ | $1$ | $ C_{1} $ | $C_{3}\times{}C_{18}$  | $\langle4,-6,1\rangle$|$\left(\begin{smallmatrix}5&-1\\4&-1\end{smallmatrix}\right)^{6}$|$\text{Y}$|$13$|
| | | $4$ | $1$ | $ C_{1} $ | $C_{6}\times{}C_{18}$  | $\langle5,-10,1\rangle$|$\left(\begin{smallmatrix}19&-2\\10&-1\end{smallmatrix}\right)^{3}$||$7$|
| | | $8$ | $2$ | $C_{2}$  | $C_{6}\times{}C_{18}$  | $\langle1,-18,1\rangle$|$\left(\begin{smallmatrix}18&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle5,-20,4\rangle$|$\left(\begin{smallmatrix}19&-4\\5&-1\end{smallmatrix}\right)^{3}$||$7$|
| $20$ | $357$ | $1$ | $2$ | $C_{2}$  | $C_{2}^{3}\times{}C_{24}$  | $\langle1,-19,1\rangle$|$\left(\begin{smallmatrix}19&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle7,-21,3\rangle$|$\left(\begin{smallmatrix}20&-3\\7&-1\end{smallmatrix}\right)^{3}$||$7$|
| $21$ | $44$ | $1$ | $1$ | $ C_{1} $ |  $C_{2}^{2}\times{}C_{24}$  | $\langle5,-8,1\rangle$|$\left(\begin{smallmatrix}22&-3\\15&-2\end{smallmatrix}\right)^{3}$||$10$|
| | | $3$ | $4$ | $C_{4}$  | $C_{2}\times{}C_{6}^{2}$  | $\langle1,-20,1\rangle$|$\left(\begin{smallmatrix}20&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle5,-24,9\rangle$|$\left(\begin{smallmatrix}22&-9\\5&-2\end{smallmatrix}\right)^{3}$||$10$|
| | | | | |  | $\langle11,-22,2\rangle$|$\left(\begin{smallmatrix}21&-2\\11&-1\end{smallmatrix}\right)^{3}$||$7$|
| | | | | |  | $\langle9,-24,5\rangle$|$\left(\begin{smallmatrix}22&-5\\9&-2\end{smallmatrix}\right)^{3}$||$10$|
| $22$ | $437$ | $1$ | $1$ | $ C_{1} $ | $C_{2}\times{}C_{120}$  | $\langle1,-21,1\rangle$|$\left(\begin{smallmatrix}21&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| $23$ | $120$ | $1$ | $2$ | $C_{2}$  | $C_{176}$  | $\langle6,-12,1\rangle$|$\left(\begin{smallmatrix}23&-2\\12&-1\end{smallmatrix}\right)^{3}$||$7$|
| | | | | |  | $\langle3,-12,2\rangle$|$\left(\begin{smallmatrix}23&-4\\6&-1\end{smallmatrix}\right)^{3}$||$7$|
| | | $2$ | $4$ | $C_{2}^{2}$  | $C_{176}$  | $\langle1,-22,1\rangle$|$\left(\begin{smallmatrix}22&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle19,-10,-5\rangle$|$\left(\begin{smallmatrix}16&5\\19&6\end{smallmatrix}\right)^{3}$||$13$|
| | | | | |  | $\langle8,-24,3\rangle$|$\left(\begin{smallmatrix}23&-3\\8&-1\end{smallmatrix}\right)^{3}$||$7$|
| | | | | |  | $\langle15,0,-8\rangle$|$\left(\begin{smallmatrix}11&8\\15&11\end{smallmatrix}\right)^{3}$||$10$|
| $24$ | $21$ | $1$ | $1$ | $ C_{1} $ | $C_{2}\times{}C_{4}\times{}C_{12}$  | $\langle1,-5,1\rangle$|$\left(\begin{smallmatrix}5&-1\\1&0\end{smallmatrix}\right)^{6}$||$7$|
| | | $5$ | $2$ | $C_{2}$  | $C_{2}^{2}\times{}C_{4}\times{}C_{12}$  | $\langle1,-23,1\rangle$|$\left(\begin{smallmatrix}23&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle17,-27,3\rangle$|$\left(\begin{smallmatrix}25&-3\\17&-2\end{smallmatrix}\right)^{3}$||$10$|
| $25$ | $572$ | $1$ | $2$ | $C_{2}$  | $C_{5}\times{}C_{40}$  | $\langle1,-24,1\rangle$|$\left(\begin{smallmatrix}24&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle13,-26,2\rangle$|$\left(\begin{smallmatrix}25&-2\\13&-1\end{smallmatrix}\right)^{3}$||$7$|
| $26$ | $69$ | $1$ | $1$ | $ C_{1} $ | $C_{2}\times{}C_{12}^{2}$  | $\langle3,-9,1\rangle$|$\left(\begin{smallmatrix}26&-3\\9&-1\end{smallmatrix}\right)^{3}$||$7$|
| | | $3$ | $3$ | $C_{3}$  | $C_{2}\times{}C_{12}^{2}$  | $\langle1,-25,1\rangle$|$\left(\begin{smallmatrix}25&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle17,-3,-9\rangle$|$\left(\begin{smallmatrix}14&9\\17&11\end{smallmatrix}\right)^{3}$||$10$|
| | | | | |  | $\langle5,-29,11\rangle$|$\left(\begin{smallmatrix}27&-11\\5&-2\end{smallmatrix}\right)^{3}$||$10$|
| $27$ | $168$ | $1$ | $2$ | $C_{2}$  | $C_{9}\times{}C_{18}$  | $\langle7,-14,1\rangle$|$\left(\begin{smallmatrix}27&-2\\14&-1\end{smallmatrix}\right)^{3}$||$7$|
| | | | | |  | $\langle11,-16,2\rangle$|$\left(\begin{smallmatrix}29&-4\\22&-3\end{smallmatrix}\right)^{3}$||$13$|
| | | $2$ | $4$ | $C_{2}^{2}$  | $C_{9}\times{}C_{18}$  | $\langle1,-26,1\rangle$|$\left(\begin{smallmatrix}26&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle7,14,-17\rangle$|$\left(\begin{smallmatrix}6&17\\7&20\end{smallmatrix}\right)^{3}$||$7$|
| | | | | |  | $\langle19,-8,-8\rangle$|$\left(\begin{smallmatrix}17&8\\19&9\end{smallmatrix}\right)^{3}$||$10$|
| | | | | |  | $\langle11,-32,8\rangle$|$\left(\begin{smallmatrix}29&-8\\11&-3\end{smallmatrix}\right)^{3}$||$10$|
| $28$ | $29$ | $1$ | $1$ | $ C_{1} $ | $C_{2}^{2}\times{}C_{6}^{2}$  | $\langle5,-7,1\rangle$|$\left(\begin{smallmatrix}6&-1\\5&-1\end{smallmatrix}\right)^{6}$|$\text{Y}$|$16$|
| | | $5$ | $2$ | $C_{2}$  | $C_{2}^{3}\times{}C_{6}^{2}$  | $\langle1,-27,1\rangle$|$\left(\begin{smallmatrix}27&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle13,-33,7\rangle$|$\left(\begin{smallmatrix}30&-7\\13&-3\end{smallmatrix}\right)^{3}$||$13$|
| $29$ | $780$ | $1$ | $4$ | $C_{2}^{2}$  | $C_{280}$  | $\langle1,-28,1\rangle$|$\left(\begin{smallmatrix}28&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle15,-30,2\rangle$|$\left(\begin{smallmatrix}29&-2\\15&-1\end{smallmatrix}\right)^{3}$||$7$|
| | | | | |  | $\langle10,-30,3\rangle$|$\left(\begin{smallmatrix}29&-3\\10&-1\end{smallmatrix}\right)^{3}$||$7$|
| | | | | |  | $\langle6,-30,5\rangle$|$\left(\begin{smallmatrix}29&-5\\6&-1\end{smallmatrix}\right)^{3}$||$7$|
| $30$ | $93$ | $1$ | $1$ | $ C_{1} $ |  $C_{2}\times{}C_{6}\times{}C_{24}$  | $\langle7,-11,1\rangle$|$\left(\begin{smallmatrix}31&-3\\21&-2\end{smallmatrix}\right)^{3}$||$10$|
| | | $3$ | $3$ | $C_{3}$  | $C_{2}\times{}C_{6}\times{}C_{24}$  | $\langle1,-29,1\rangle$|$\left(\begin{smallmatrix}29&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle19,-1,-11\rangle$|$\left(\begin{smallmatrix}15&11\\19&14\end{smallmatrix}\right)^{3}$||$10$|
| | | | | |  | $\langle7,-33,9\rangle$|$\left(\begin{smallmatrix}31&-9\\7&-2\end{smallmatrix}\right)^{3}$||$10$|
| $31$ | $56$ | $1$ | $1$ | $ C_{1} $ | $C_{10}\times{}C_{30}$  | $\langle2,-8,1\rangle$|$\left(\begin{smallmatrix}31&-4\\8&-1\end{smallmatrix}\right)^{3}$||$7$|
| | | $2$ | $2$ | $C_{2}$  | $C_{10}\times{}C_{30}$  | $\langle8,-16,1\rangle$|$\left(\begin{smallmatrix}31&-2\\16&-1\end{smallmatrix}\right)^{3}$||$7$|
| | | | | |  | $\langle11,2,-5\rangle$|$\left(\begin{smallmatrix}13&10\\22&17\end{smallmatrix}\right)^{3}$||$13$|
| | | $4$ | $4$ | $C_{4}$  | $C_{10}\times{}C_{30}$  | $\langle1,-30,1\rangle$|$\left(\begin{smallmatrix}30&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle13,-34,5\rangle$|$\left(\begin{smallmatrix}32&-5\\13&-2\end{smallmatrix}\right)^{3}$||$10$|
| | | | | |  | $\langle25,-36,4\rangle$|$\left(\begin{smallmatrix}33&-4\\25&-3\end{smallmatrix}\right)^{3}$||$13$|
| | | | | |  | $\langle5,-34,13\rangle$|$\left(\begin{smallmatrix}32&-13\\5&-2\end{smallmatrix}\right)^{3}$||$10$|
| $32$ | $957$ | $1$ | $2$ | $C_{2}$  | $C_{2}\times{}C_{16}^{2}$  | $\langle1,-31,1\rangle$|$\left(\begin{smallmatrix}31&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle11,-33,3\rangle$|$\left(\begin{smallmatrix}32&-3\\11&-1\end{smallmatrix}\right)^{3}$||$7$|
| $33$ | $1020$ | $1$ | $4$ | $C_{2}^{2}$  | $C_{2}\times{}C_{120}$  | $\langle1,-32,1\rangle$|$\left(\begin{smallmatrix}32&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle17,-34,2\rangle$|$\left(\begin{smallmatrix}33&-2\\17&-1\end{smallmatrix}\right)^{3}$||$7$|
| | | | | |  | $\langle23,-36,3\rangle$|$\left(\begin{smallmatrix}34&-3\\23&-2\end{smallmatrix}\right)^{3}$||$10$|
| | | | | |  | $\langle29,-40,5\rangle$|$\left(\begin{smallmatrix}36&-5\\29&-4\end{smallmatrix}\right)^{3}$||$16$|
| $34$ | $1085$ | $1$ | $2$ | $C_{2}$  | $C_{2}\times{}C_{288}$  | $\langle1,-33,1\rangle$|$\left(\begin{smallmatrix}33&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle7,-35,5\rangle$|$\left(\begin{smallmatrix}34&-5\\7&-1\end{smallmatrix}\right)^{3}$||$7$|
| $35$ | $8$ | $1$ | $1$ | $ C_{1} $ | $C_{6}\times{}C_{12}$  | $\langle2,-4,1\rangle$|$\left(\begin{smallmatrix}3&-1\\2&-1\end{smallmatrix}\right)^{12}$|$\text{Y}$|$13$|
| | | $2$ | $1$ | $ C_{1} $ | $C_{2}\times{}C_{6}\times{}C_{12}$  | $\langle1,-6,1\rangle$|$\left(\begin{smallmatrix}6&-1\\1&0\end{smallmatrix}\right)^{6}$||$7$|
| | | $3$ | $1$ | $ C_{1} $ | $C_{2}\times{}C_{6}\times{}C_{24}$  | $\langle7,-10,1\rangle$|$\left(\begin{smallmatrix}37&-4\\28&-3\end{smallmatrix}\right)^{3}$||$13$|
| | | $4$ | $1$ | $ C_{1} $ | $C_{2}\times{}C_{6}\times{}C_{24}$  | $\langle4,-12,1\rangle$|$\left(\begin{smallmatrix}35&-3\\12&-1\end{smallmatrix}\right)^{3}$||$7$|
| | | $6$ | $2$ | $C_{2}$  | $C_{2}\times{}C_{6}\times{}C_{24}$  | $\langle9,-18,1\rangle$|$\left(\begin{smallmatrix}35&-2\\18&-1\end{smallmatrix}\right)^{3}$||$7$|
| | | | | |  | $\langle4,-20,7\rangle$|$\left(\begin{smallmatrix}37&-14\\8&-3\end{smallmatrix}\right)^{3}$||$10$|
| | | $12$ | $4$ | $C_{4}$  | $C_{2}\times{}C_{6}\times{}C_{24}$  | $\langle1,-34,1\rangle$|$\left(\begin{smallmatrix}34&-1\\1&0\end{smallmatrix}\right)^{3}$||$4$|
| | | | | |  | $\langle16,8,-17\rangle$|$\left(\begin{smallmatrix}13&17\\16&21\end{smallmatrix}\right)^{3}$||$13$|
| | | | | |  | $\langle4,-36,9\rangle$|$\left(\begin{smallmatrix}35&-9\\4&-1\end{smallmatrix}\right)^{3}$||$7$|
| | | | | |  | $\langle28,-12,-9\rangle$|$\left(\begin{smallmatrix}23&9\\28&11\end{smallmatrix}\right)^{3}$||$13$|
