Download Link: https://assignmentchef.com/product/solved-ecse343-assignment-2
<br>



Please submit this <strong>.mlx</strong> file along with the <strong>PDF</strong> copy of this file.

<strong>Question 1: </strong>

<ol>

 <li>Prove that .</li>

</ol>

Ax = 0

Ax = 0, so X belongs to N(          A)  hence, N(A) belongs to N(A); also x belongs to N(A)  so Ax = 0

Ax =0

Ax = 0, so x belongs to N(A)  because N(A) belongs to N(A)

so N(A) = N(A) and thier dimensions are the sane  finally

<ol>

 <li>Given two orthonormal matrices Uand Q, show that product of these matrices, UQ, is also orthonormal.</li>

</ol>

U             =            U = I Q           =              Q = I so            (UQ) =    UQ = I

in this way, it is orthonormal

<ol>

 <li>c) If A is an invertible matrix and Q is an orthonormal matrix, show that</li>

</ol>

Hint:  (the ratio of the largest and smallest eigen values).

A is a intvertiable matrix and ! is an orthonormal matrix, in this way, det(A) n != 0 also det(A) is the product of eigenvalues we know that eigenvalues of A are non-zero

Q = I

det (           )=det(

(detQ)^2 = 1 so that det(Q) = +1/-1

det(QA) = det(Q)(detA) = (+/- 1)*(product of eigenvalues of A)

K(QA) = =K(A)

<strong>Question 2: </strong>

(a) <strong>Cholesky factorization </strong>is a popular matrix factorization algorithm analogous to LU decompostion for positive definite matrices. In this case a matrix  is factored as




Where  is a lower triangular matrix with real and positive diagonal enteries.

Note that there are a number of tests/definitions for positive definite matrices. A matrix, Mis positive definite if and only if it satisfies the following two conditions,

<ol>

 <li>The matrix Mmust be sysmterical, i.e.</li>

 <li>for  all  .</li>

</ol>

An alternative definition/test is that a matrix is positive definite if and only if it can be written as




It can be shown (similarly to how we derived the Doolittle algorithm) that the entries of matrix L are given by the following expression

The cost of Cholesky factorization algorithm is rougly half than the LU decompostion. Your task is to implement the Cholesky decomposition alogorithm in the cell below. Note that for Cholesky decomposition pivoting is not required

The above equation consists of the triangular systems, the solution <strong>x </strong>can be obtained by  first solving  using the forward substitution, then by solving the   using backward subsitution. Implement the Cholesky solver for least squares in the function below.

You can use your Forward/Backward Substitution from Assignment 1 for bonus marks. Include your code for forward/backward substitution in the Appendix.

% in the end

<ol>

 <li>d) To test the Cholesky Least Squares solver we numerical algorithm we will use a square matrix. Run the cell below to compare the solution.</li>

</ol>

<strong>Question 3: </strong>

Polynomial fitting is one of the many applications of Least-Squares Approximation. Given the experimental data, shown in Table below,

we can fit a degree polynomial,   . To find the coefficicients,  , we can solve the linear system shown below.

Where M is a Vandermonde matrix (known to be ill-conditioned). The system                 can then be solved to compute the coefficients in vector

<ol>

 <li>Your first task is to write a function that takes the input data vector x and degree of the polynomial, then returns the Vadermonde</li>

 <li>The linear system can be solved using different approaches. In the code below, use the Cholesky decompositon and QR decomposition function to find the polynomial coefficients for the 1st degree polynomial. Explain the results obtained in figure 1.</li>

</ol>

% Load the Input Data from the provided .mat file

<ol start="2">

 <li>Now use the cell below to use the Cholesky decompositon and QR decomposition function to find the polynomial coefficients for the 14th degree polynomial. <strong>Explain the results obtained in figure 2.</strong></li>

 <li>In this part we will use limited amount of data, we will use only<strong> 15 data points</strong>. Use Cholesky decompositon and QR decomposition function to find the polynomial coefficients for the 14th degree polynomial. <strong>Comment on the results obtained in figure 3, which algorithm performs better and why?</strong></li>

</ol>

<strong>Question 4:</strong>

<strong> </strong>In this question we will develop the different algorithms to implment QR factorization.

<ol>

 <li>Use the cell below to implement a function that computes Gram-Schmidt based QR decompositon of a input matrix A.</li>

 <li>Use the cell below to implement a function that uses the Modiefied Gram-Schmidt approach to compute the QR decompositon of  the input matrix A.</li>

 <li>Use the cell below to implement a function that uses the Householder approach to compute the QR decompositon of  the input matrix A.</li>

 <li>d) We know that Q matrix obtained using QR decompostion should have orthonormal columns i.e.</li>

</ol>

numerical errors, In this section we will test the functions written above to compare the error between the above implemeneted approaches. In oder to do this we will use a  HIlbert Matrix, the entries of Hibert Matrix are given by

In the cell below we compute the QR decomposition of Hilbert matrices of size ranging from 2 to 16, then we compute the norm of error  to assess the performance of the lorithms written in part (a), ( b), and (c). Comment on the results obtained. for n=2:16

a= hilb(n); % n is the size of the matrix

[q_gs,r1] = gschmidt(a); %using algorithm classical Gram–Schmidt [q_mgs,r2] = mgschmidt(a); %using algorithm modified Gram-Schmidt

[q_hh,r3] = householder(a); % using house holder err_gs(n-1) = norm(q_gs’*q_gs – eye(n),’fro’); err_mgs(n-1) = norm(q_mgs’*q_mgs – eye(n),’fro’); err_hh(n-1) = norm(q_hh’*q_hh – eye(n),’fro’); end

1 0000    0 0510    0 3333




0.2000    0.0479    0.0006       0.2500    0.0515   -0.0016   -0.0005    0.0000    0.0000   -0.0000    0.0000   -0.0000    0.0.0005   0.0000    0.0000   0.0000   0.0000    0.0000   0.     0.1667    0.0438    0.0019   -0.0003   -0.0000    0.0000    0.0000   -0.0000   -0.0000    0.    0.2000    0.0479    0.0006   -0.0005   -0.0000    0.0000   -0.0000   -0.0000    0.0000   -0.     0.1429    0.0401    0.0027   -0.0002   -0.0000   -0.0000    0.0000    0.0000   -0.0000    0.    0.1667    0.0438    0.0019   -0.0003   -0.0000    0.0000    0.0000   -0.0000   -0.0000    0.     0.1250    0.0369    0.0031   -0.0001   -0.0000   -0.0000    0.0000    0.0000   -0.0000   -0.    0.1429    0.0401    0.0027   -0.0002   -0.0000   -0.0000    0.0000    0.0000   -0.0000    0.     0.1111    0.0340    0.0034    0.0001   -0.0000   -0.0000    0.0000    0.0000    0.0000   -0.    0.1250    0.0369    0.0031   -0.0001   -0.0000   -0.0000    0.0000    0.0000   -0.0000   -0.     0.1000    0.0315    0.0035    0.0001   -0.0000   -0.0000   -0.0000    0.0000    0.0000    0.    0.1111    0.0340    0.0034    0.0001   -0.0000   -0.0000    0.0000    0.0000    0.0000   -0.

0.1000    0.0315    0.0035    0.0001   -0.0000   -0.0000   -0.0000    0.0000    0.0000    0.

<ol start="3">

 <li>e) Finally, now we have implemented various different versions of QR decompositon. Use above versions of the the QR transform and use this to fit the polynomial in question 3. We want to fit the polynomial of degree 14 through the data with 15 points. Clearly show all the plots (with legends) and comment on the results obtained.</li>

</ol>

load(‘Polynomial_Fitting_Data.mat’);

v =

0.0000   -0.0000         0         0         0         0         0         0         0

0.0000   -0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.     0.0000   -0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.     0.0000   -0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.     0.0000   -0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.     0.0000   -0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0. v =   0.0000   -0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.

0.0000   -0.0000   -0.0000         0         0         0         0         0         0                                                  0.000     0.000     0.000     0.000     0.000      0.0001   0.

0.0000   -0.0000   -0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.              1      4     0.0000   -0.0000   -0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.           1      2      1     0.0000   -0.0000   -0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.     0.0000   -0.0000   -0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.     0.0000   -0.0000   -0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0. v =   0.0000   -0.0000   -0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.     0.0000         0   -0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0001    0.     0.0000    0.0000   -0.0000    0.0000    0.0000    0.0000    0.0000    0.0001    0.0004    0.

0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0001    0.0002    0.0010    0.

78/92

0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.     0.0000   -0.0000    0.0000   -0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0000    0.     0.0000   -0.0000    0.0000    0.0000   -0.0000    0.0000   -0.0000   -0.0000    0.0000   -0.

0.0000   -0.0000   -0.0000    0.0000   -0.0000    0.0000    0.0000   -0.0000    0.0000    0.

v =     0.0000   -0.0000   -0.0000    0.0000   -0.0000   -0.0000    0.0000   -0.0000   -0.0000    0.

0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.    0.0000   -0.0000   -0.0000    0.0000    0.0000   -0.0000    0.0000    0.0000   -0.0000   -0.     0.0000   -0.0000    0.0000   -0.0000   -0.0000    0.0000   -0.0000    0.0001   -0.0001    0.    0.0000   -0.0000   -0.0000    0.0000    0.0000   -0.0000   -0.0000    0.0000    0.0000   -0.     0.0000   -0.0000    0.0000    0.0000   -0.0000    0.0000   -0.0000   -0.0000    0.0001   -0.    0.0000         0   -0.0000    0.0000    0.0000    0.0000   -0.0000    0.0000    0.0000    0.     0.0000   -0.0000   -0.0000    0.0000   -0.0000    0.0000    0.0000   -0.0001    0.0000    0.    0.0000    0.0000   -0.0000   -0.0000    0.0000    0.0000   -0.0000   -0.0000    0.0000    0.     0.0000   -0.0000   -0.0000    0.0000   -0.0000   -0.0000    0.0000   -0.0000   -0.0001    0.    0.0000    0.0000   -0.0000   -0.0000    0.0000    0.0000    0.0000   -0.0000   -0.0000    0.     0.0000   -0.0000   -0.0000    0.0000    0.0000   -0.0000    0.0000    0.0000   -0.0001   -0. v =   0.0000   -0.0000   -0.0000    0.0000    0.0000   -0.0000   -0.0000    0.0000    0.0000   -0.

0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.                               –                                         –                         1

0.0000   -0.0000    0.0000   -0.0000   -0.0000    0.0000   -0.0000    0.0001   -0.0001    0. – – 0 0     0.0000   -0.0000    0.0000    0.0000   -0.0000    0.0000   -0.0000   -0.0000    0.0001   -0. – –  –     0.0000   -0.0000   -0.0000    0.0000   -0.0000    0.0000    0.0000   -0.0001    0.0000    0.     0.0000   -0.0000   -0.0000    0.0000   -0.0000   -0.0000    0.0000   -0.0000   -0.0001    0.     0.0000   -0.0000   -0.0000    0.0000    0.0000   -0.0000    0.0000    0.0000   -0.0001   -0. v =   0.0000   -0.0000   -0.0000    0.0000    0.0000   -0.0000   -0.0000    0.0000    0.0000   -0.

0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.0000    0.0000   -0.                               –                                         –                         1

0.0000   -0.0000    0.0000   -0.0000   -0.0000    0.0000   -0.0000    0.0001   -0.0001    0. – – 0 0     0.0000   -0.0000    0.0000    0.0000   -0.0000    0.0000   -0.0000   -0.0000    0.0001   -0. – –  –     0.0000   -0.0000   -0.0000    0.0000   -0.0000    0.0000    0.0000   -0.0001    0.0000    0.     0.0000   -0.0000   -0.0000    0.0000   -0.0000   -0.0000    0.0000   -0.0000   -0.0001    0.     0.0000   -0.0000   -0.0000    0.0000    0.0000   -0.0000    0.0000    0.0000   -0.0001   -0. v =   0.0000   -0.0000   -0.0000    0.0000    0.0000   -0.0000   -0.0000    0.0000    0.0000   -0.

0.0000   -0.0000    0.0001   -0.0001    0.0002   -0.0002    0.0002   -0.0002    0.0002   -0.                               –     0          0          0          0 –      0          0          1

0.0000   -0.0000    0.0000   -0.0000   -0.0001    0.0002   -0.0005    0.0007   -0.0009    0. – 0 0 0 – 0 0     0.0000   -0.0000    0.0000    0.0000   -0.0001    0.0002   -0.0001   -0.0003    0.0010   -0. – – 0 0 0 0 – 01     0.0000   -0.0000   -0.0000    0.0001   -0.0001    0.0000    0.0003   -0.0006    0.0004    0.     0.0000   -0.0000   -0.0000    0.0001   -0.0000   -0.0002    0.0003   -0.0001   -0.0007    0.     0.0000   -0.0000   -0.0000    0.0001    0.0000   -0.0002    0.0001    0.0005   -0.0007   -0. v =   0.0000   -0.0000   -0.0000    0.0000    0.0001   -0.0001   -0.0002    0.0005    0.0003   -0.

0.0000   -0.0000    0.0001   -0.0001    0.0002   -0.0002    0.0002   -0.0002    0.0002   -0.                               –     0          0          1          0 –      3          0          8

0.0000   -0.0000    0.0000   -0.0000   -0.0001    0.0002   -0.0005    0.0007   -0.0009    0. – 1 2 – 5 3     0.0000   -0.0000    0.0000    0.0000   -0.0001    0.0002   -0.0001   -0.0003    0.0010   -0. – – 1 0 5 – 07     0.0000   -0.0000   -0.0000    0.0001   -0.0001    0.0000    0.0003   -0.0006    0.0004    0.     0.0000   -0.0000   -0.0000    0.0001   -0.0000   -0.0002    0.0003   -0.0001   -0.0007    0.     0.0000   -0.0000   -0.0000    0.0001    0.0000   -0.0002    0.0001    0.0005   -0.0007   -0. v =   0.0000   -0.0000   -0.0000    0.0000    0.0001   -0.0001   -0.0002    0.0005    0.0003   -0.

0.0001   -0.0003    0.0007   -0.0011    0.0016   -0.0021    0.0025   -0.0025    0.0022   -0.0      0 –                            0        00         01          00 –     03        00         08

0.0001   -0.0003    0.0004   -0.0002   -0.0007    0.0025   -0.0049    0.0075   -0.0090    0.0 0 – 0 0 1 01 02 – 0 03     0.0001   -0.0002    0.0001    0.0004   -0.0014    0.0021   -0.0009   -0.0033    0.0099   -0.0 0 – 0 – 1 00 02 1 05 – 07     0.0001   -0.0002   -0.0001    0.0007   -0.0012    0.0001    0.0030   -0.0059    0.0037    0.     0.0001   -0.0001   -0.0002    0.0008   -0.0004   -0.0016    0.0034   -0.0006   -0.0074    0.     0.0001   -0.0001   -0.0003    0.0006    0.0004   -0.0021    0.0009    0.0048   -0.0066   -0. v =   0.0001   -0.0000   -0.0004    0.0003    0.0010   -0.0014   -0.0022    0.0048    0.0030   -0.

1.0000   -3.2776    6.6502  -11.2097   16.4958  -21.4550   24.6579  -24.8673   21.7346  -16.0  1      0 -0 00 4  0 0 00 0 0012  0 000  -0 0034  0 0000     0 0084  0     1.0000   -2.8094    3.8001   -1.6014   -7.0696   24.5200  -49.3159   74.6019  -90.0434   87.0   1 0 0 00 -0 0 4 0 0 03 0 010 0 0014 -0 0022 -0 0 48  0 0 0     1.0000   -2.3411    1.3885    4.3114  -14.3206   20.9835   -9.4838  -32.5188   99.1194 -160.0      1 0 000       -0 0003 -0 0006  0 00 4 0021 0 0009 -0 004 -0 0066

1.0000   -1.8729   -0.5846    7.1446  -11.6015    0.9431   30.3482  -59.2989   37.4982   61.

1.0000   -1.4047   -2.1193    7.5142   -4.1034  -16.0966   33.9693   -5.7386  -74.2799  117.

1.0000   -0.9365   -3.2155    6.0360    4.1363  -21.4336    8.6217   47.8217  -65.6815  -44.4

A3 = (q_gs*r2)y;

Warning: Matrix is close to singular or badly scaled. Results may be inaccurate. RCOND =  6.045334e-25.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_new = 0:7/15:7;    x_new = x_new’;




fa = zeros(length(x_new),1);    fb = zeros(length(x_new),1);    fc = zeros(length(x_new),1);







L(i,j)=M(i,j)/L(j,j);     elseif i&gt;j

L(i,j)=1/L(j,j)*M(i,j);     elseif i&lt;j       L(i,j)=0;     end   end end end

% question 2 b function x = LeastSquareSolver_Cholesky(A,b)

M = A’*A;

L = CholeskyDecomposition(M); y = forward_sub(L,b);   x = backward_sub(L,b); end function y = forward_sub(L,b) l = length(b); y = zeros(l,1); y(1,1) = b(1)./L(1,1); for j = 2:l

y(j,1) = (b(j)-sum(L(j,1:j-1)*y(1:j-1,1)))./L(j,j);

end end

function x = backward_sub(L,b)

l = length(b); x = zeros(1,l); x(1,l)=b(end)./L(l,l); for i = l-1:-1:1

temp = 1/L(i,i).*(b(i)-sum(L(i,i+1:end).*x(i+1:end)));

x(1,i) = temp;

end x=x’ end

%question 3a

function M = PolynomialMatrix(x,n)

powers = 0:n; for i = 1:length(powers)     M(:,i) = x.^powers(i);

end end

%question 4a

function [Q,R] = gschmidt(A)

[n,p] = size(A);

Q = zeros(n,p);     R = zeros(p,p);     for i = 1:p         Q(:,i) = A(:,i);         if i ~= 1

R(1:i-1,i) = Q(:,i-1)’*Q(:,i);

Q(:,i) = Q(:,i) – Q(:,1:i-1)*R(1:i-1,i);

end

R(i,i) = norm(Q(:,i));

Q(:,i) = Q(:,i)/R(i,i);      end end   %question 4b

function [Q,R] = mgschmidt(A)

[n,p] = size(A);

<ul>

 <li>= zeros(n,p);</li>

 <li>= zeros(p); v = zeros(p);     for i = 1:n         v(i,:) = A(i,:);     end         for i = 1:n             R(i,i) = normest(v(:,i));</li>

</ul>

Q(:,i) = v(:,i)/R(i,i);        for j = (i+1):p