---

---

# Курсовая работа

## 1. **Условие задачи**

### Вариант 43

МКЭ для двумерной краевой задачи для эллиптического уравнения в декартовой системе
координат. Базисные функции биквадратичные на прямоугольниках. Краевые условия
всех типов. Коэффициент диффузии  разложить по билинейным базисным функциям.
Матрицу СЛАУ генерировать в разреженном строчном формате. Для решения СЛАУ
использовать МСГ или ЛОС с неполной факторизацией.

## 2. **Постановка задачи**

Эллиптическая краевая задача для функции $u$ определяется дифференциальным уравнением:

$$
-div(\lambda \space grad \space u) + \gamma u = f, \qquad (2.1)
$$

заданным в некоторой области $\Omega$ с границей $S = S_1 \cup S_2 \cup S_3$, и краевыми условиями:

$$
u\bigg|_{S_1} = u_g,
$$

$$
\lambda \space \frac{\partial u}{\partial n} \bigg|_{S_2} = \theta,
$$

$$
\lambda \space \frac{\partial u}{\partial n} \bigg|_{S_3} + \beta (u\big|_{S_3} - u_{\beta}) = 0.
$$

## 3. **Теоретическая часть**

### 3.1 **Вариационная постановка в форме уравнения Галёркина**

Запишем для нашей краевой задачи эквивалентную вариационную постановку в форме уравнения Галёркина

$$
R(u) = 0,
$$

$$
R(u) = -div(\lambda \space grad \space u) + \gamma u - f
$$

где $R(u)$ - невязка уравнения (2.1).

Потребуем, чтобы невязка $R(u)$ была ортогональна некоторому пространству $\Phi$ функций $v$, которое будем называть _пространством пробных функций_. Ортогональность в том смысле, что скалярное произведение $(R(u),v)$ в пространстве $L_2(\Omega)$ равно 0, т.е.

$$
\int\limits_{\Omega} (-div(\lambda \space grad \space u) + \gamma u - f)vd\Omega = 0,\quad\forall v \in \Phi \qquad (3.1)
$$

Воспользуемся _формулой Грина_ для преобразования интеграла $\int\limits_{\Omega} (-div(\lambda \space grad \space u)vd\Omega$ в (3.1). Получим:

$$
{\int\limits_{\Omega} \lambda \space grad \space u \cdot grad \space v \space d\Omega}\space - {\int\limits_{S}\lambda \space {\frac{\partial u}{\partial n}} vdS } \space + {\int\limits_{\Omega} (\gamma u - f )vd\Omega} = 0,\quad\forall v \in \Phi \qquad (3.2)
$$

По условию граница $S = S_1 \cup S_2 \cup S_3$, на которых заданы краевые условия трех типов. Преобразуем интеграл ${ \int\limits_{S}\lambda \space {\frac{\partial u}{\partial n}} vdS }$ в соответствие с краевыми условиями:

1. ${ \int\limits_{S_1}\lambda \space {\frac{\partial u}{\partial n}} vdS } = 0$.
   Так как на границе $S_1$ задается только значение функции $u$, то нужно потребовать, чтобы пространство пробных функций $\Phi$ содержало только те функции, которые принимают нулевое значение на границе $S_1$.
   Далее такие функции будем называть $v_0$.

2. ${ \int\limits_{S_2}\lambda \space {\frac{\partial u}{\partial n}} vdS } = {\int\limits_{S_2} \theta v dS}$

3. $ { \int\limits_{S_3}\lambda \space {\frac{\partial u}{\partial n}} vdS } = -{ \int\limits_{S_3} \beta(u-u_{\beta})vdS }$

Таким образом, получаем вариационное уравнение вида:

$$
{\int\limits_{\Omega} \lambda \space grad \space u \cdot grad \space v_0 d \Omega} \space + {\int\limits_{\Omega} \gamma u v_0 d \Omega} \space + {\int\limits_{S_3} \beta u v_0 dS } =

\\

= {\int\limits_{\Omega} f v_0 d \Omega} \space + {\int\limits_{S_2} \theta v_0 dS} \space + {\int\limits_{S_3} \beta u_{\beta} v_0 dS}, \quad \forall v_0 \in H_0^1 \qquad \qquad (3.3)
$$

где $H_0^1$ - пространство функций, имеющих суммируемые с квадратом производные и равных нулю на границе $S_1$

При этом будем считать, что $u \in H_g^1$, где $H_g^1$ - множество функций, имеющих суммируемые с квадратом первые производные и удовлетворяющих только первым краевым условиям на границе $S_1$ (по аналогии с пространством $H_0^1$).

### 3.2 **Дискретизация**

Для построения конечноэлементных аппроксимаций по методу Галёркина пространства $H_g^1$ и $H_0^1$ заменяются конечномерными пространствами $V_g^h$ и $V_0^h$.

В МКЭ базисом пространства $V^h$ является набор финитных кусочно-полиномиальных функций $\psi_i, i \in \overline{1,n}$.

> Финитные функции - это такие функции, каждая из которых отлична от нуля лишь на нескольких малых подобластях $\Omega_i$ расчетной области $\Omega$.

Получим аппроксимацию уравнения (3.3) заменой функции $u$ аппроксимирующей ее функцией $u^h$, а функцию $v_0$ заменим функцией $v_0^h$

Любая функция $v_0^h \in V_0^h$ может быть представлена в виде линейной комбинации:

$$
v_0^h = \sum\limits_{i \in N_0} q_i^v \psi_i, \qquad (3.4)
$$

Преобразуем уравнение (3.3) исходя из уравнения (3.4):

$$
{\int\limits_{\Omega} \lambda \space grad \space u^h \cdot grad \space \psi_i d \Omega} + {\int\limits_{\Omega} \gamma u^h \psi_i d \Omega} + {\int\limits_{S} \beta u^h \psi_i dS} = 
\\

= {\int\limits_{\Omega} f \psi_i d \Omega} + {\int\limits_{S_2} \theta \psi_i dS} + {\int\limits_{S_3} \beta u_{\beta} \psi_i dS}, \quad i \in N_0 \qquad \qquad \qquad (3.5)
$$

Таким образом, МКЭ-решение $u^h$ удовлетворяет системе уравнений (3.5).

В свою очередь, функция $u^h$ так же может быть представлена в виде линейной комбинации базисных функций пространства $V^h$:

$$
u^h = \sum\limits_{j = 1} ^ {n} q_i \psi_i, \qquad (3.6)
$$

Подставляя (3.6) в (3.5), получаем итоговое конечноэлементное СЛАУ _для компонент_ $q_j$ вектора весов $q = (q_1,...,q_n)^T$ *с индексами* $j \in N_0$:

$$
\sum\limits_{j = 1}^n {\Biggl ( } { {\int\limits_{\Omega} \lambda \space grad \space \psi_j \cdot grad \space \psi_i d \Omega} + {\int\limits_{\Omega} \gamma \psi_j \psi_i d \Omega} + {\int\limits_{S_3} \beta u_{\beta} \psi_i dS} }  {\Biggl )} \space q_j =
\\
= {\int\limits_{\Omega} f \psi_i d \Omega} + {\int\limits_{S_2} \theta \psi_i dS} + {\int\limits_{S_3} \beta u_{\beta} \psi_i dS}, \quad i \in N_0, \qquad (3.7)
$$

### 3.3 **Решение задачи на прямоугольной сетке**

По условию базисные функции являются ***биквадратичными***, то есть функция принимает вид:

$$
\psi(x,y) = a_0 + a_1 x + a_2 y + a_3 xy + a_4 x^2 + a_5 y^2 +

\\

+ \space a_6 x^2 y + a_7 xy^2 + a_8 x^2 y^2 
$$

Для определения явного вида биквадратичных базисных функций на прямоугольном конечном элементе представим каждую базисную функцию в виде произведения двух одномерных квадратичных базисных функций координат $x$ и $y$, то есть:

$$
\psi(x,y) = X(x)Y(y) \qquad (3.8)
$$

Квадратичный полином задается тремя точками. Возьмем нормированный конечный элемент $[0;1] \times [0;1]$ и точки:

$$
x_1 = 0 \quad x_2 = 1/2 \quad x_3 = 1
$$

$$
y_1 = 0 \quad y_2 = 1/2 \quad y_3 = 1
$$

Таким образом образом образуется _9 узлов_.
Каждому узлу поставим в соответствие квадратичную базисную функцию:

$$
X_1 = 2(x - 1/2)(x - 1)
$$

$$
X_2 = -4 x (x - 1)
$$

$$
X_3 = 2 x (x - 1/2)
$$

Функции $Y(y)$ строятся аналогично.

Тогда локальные базисные функции на данном нормированном конечном элементе имеют вид:

$$
\hat{\psi_1}(x,y) = X_1(x) Y_1(y)

\\

\hat{\psi_2}(x,y) = X_2(x) Y_1(y)

\\

\hat{\psi_3}(x,y) = X_3(x) Y_1(y)
$$

$$
\hat{\psi_4}(x,y) = X_1(x) Y_2(y)

\\

\hat{\psi_5}(x,y) = X_2(x) Y_2(y)

\\

\hat{\psi_6}(x,y) = X_3(x) Y_2(y)
$$

$$
\hat{\psi_7}(x,y) = X_1(x) Y_3(y)

\\

\hat{\psi_8}(x,y) = X_2(x) Y_3(y)

\\

\hat{\psi_9}(x,y) = X_3(x) Y_3(y)
$$

Коэффициент диффузии $\lambda$ разложим по билинейным базисным функциям, т.е. :

$$
\lambda = \sum\limits_{i = 1}^4 \phi_i \lambda_i \qquad (3.9)
$$

$$
\phi_1^{\lambda} = xy;
$$

$$
\phi_2^{\lambda} = (1-x)y;
$$

$$
\phi_3^{\lambda} = x(1-y);
$$

$$
\phi_4^{\lambda} = (1-x)(1-y);
$$

### 3.4 **Построение матриц жесткости и масс**

СЛАУ задана в следующем виде:

$$
Aq = b
$$

и ее компоненты вычисляются следующим образом:

$$
A_{ij} = {\int\limits_{\Omega} \lambda \space grad \space \psi_j \cdot grad \space \psi_i d \Omega} \space + {\int\limits_{\Omega} \gamma \psi_j \psi_i d \Omega} \space +

\\

+ {\int\limits_{S_3} \beta \psi_j \psi_i dS }, \quad i \in N_0 \qquad (3.10)
$$

Но данный способ вычисления компонент глобальной матрицы А в многомерных задачах крайне неудобен, поэтому поступают следующим образом:

- разбивают расчетную область, т.е. $\Omega = \bigcup\limits_k \Omega_k$

- вычисляются суммы по _k_ соответствующих интегралов, т.е:

$$
G_{ij} = \sum\limits_k{\int\limits_{\Omega_k} \lambda \space grad \space \psi_j \cdot grad \space \psi_i d \Omega}, \quad i \in N_0
$$

$$
M_{ij} = \sum\limits_k{{\int\limits_{\Omega_k} \gamma \psi_j \psi_i d \Omega}}, \quad i \in N_0
$$

$$
M_{ij}^{S_3} = \sum\limits_l{\int\limits_{S_3^l} \beta \psi_j \psi_i dS }, \quad i \in N_0
$$

где $G$ - матрица жёсткости, а $M$ - матрица массы.

Локальная матрица $\hat A_{ij}$ вычисляется так:

$$
\hat A_{ij} = \hat G_{ij} + \hat M_{ij} 
$$

Локальные компоненты матриц жесткости и массы вычисляются следующим образом:

$$
\hat G_{ij} = {\int\limits_{\Omega_k} \lambda {\Biggl(} {\frac{\partial \hat\psi_i}{\partial x}}{\frac{\partial \hat\psi_j}{\partial x}} + {\frac{\partial \hat\psi_i}{\partial y}}{\frac{\partial \hat\psi_j}{\partial y}} {\Biggl)} dxdy} \qquad (3.11)
$$

$$
\hat M_{ij} = {\int\limits_{\Omega_k} \gamma \hat\psi_i \hat\psi_j dxdy}. \qquad (3.12)
$$

Локальные компоненты вектора правой части вычисляются так:

$$
\hat b_i = \int\limits_{\Omega_k} f \hat\psi_i dxdy
$$

В двумерной задаче для построения матриц жесткости и масс вычислять локальные компоненты соответствующих матриц нужно через локальные компоненты матриц одномерных квадратичных элементов. Происходит это так:

$$
\hat G_{ij} = 
{\bar \lambda {\Big (G_{\mu(i) \mu(j)}^x M_{\nu(i) \nu(j) }^y}
+ M_{\mu(i) \mu(j)}^x G_{\nu(i) \nu(j) }^y \Big )},
$$

$$
\hat M_{ij} = \bar \gamma M_{\mu(i) \mu(j)}^x M_{\nu(i) \nu(j)}^y
$$

http://www.astronet.ru/db/msg/eid/latex%20/greec.html
https://math.meta.stackexchange.com/questions/5020/mathjax-basic-tutorial-and-quick-reference
https://stackedit.io/app#







$$
kok
$$