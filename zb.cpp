#include <stdio.h>
double Lagrange(double a[], double b[], double x);
double dEt2(double a1, double b1, 
            double a2, double b2);
int main()
{
    double x[27] = { 0,0.01,0.02,0.04,0.06,0.08,0.1,0.14,0.18,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.894,0.9,0.95,1 };
    double y[27] = { 0,0.11,0.17,0.27,0.34,0.392,0.43,0.482,0.513,0.525,0.551,0.575,0.595,0.614,0.635,0.657,0.678,0.698,0.725,0.755,0.785,0.82,0.855,0.894,0.898,0.942,1 };
    int i = 0, N=1;
    double R, /*R2, x1, y1,*/ xd, xw, xf, q, xq, yq, temp;
    double a[3], b[3];
    xd = 0.8408;
    xw = 0.0408;
    xf = 0.2001;
    double yn[100] = { 0 }, xn[100] = { 0 }; yn[0] = xd, xn[0] = xd;
    R = (110 - 22) / 22.0;
    q = 1.139;
    printf("请输入R的值！");
    scanf("%lf", &R);
    for (N = 1; xn[N - 1] > 0.0381; N++)    //逐板计算法求全回流理论板数
    {
        yn[1] = 0.8662;
        for (i = 0; yn[N] > y[i]; i++);
        a[0] = y[i - 1], a[1] = y[i], a[2] = y[i + 1];
        b[0] = x[i - 1], b[1] = x[i], b[2] = x[i + 1];
        xn[N] = Lagrange(a, b, yn[N]);
        yn[N + 1] = xn[N];
    }
    printf("全回流理论板数：%d\n", N - 1);
    for (i = 1; i < 100; i++)
    {
        yn[i] = 0;
        xn[i] = 0;
    }
    //数组归零
    for (N = 1; xn[N - 1] > xf; N++)        //逐板计算法求精馏段理论板数
    {
        yn[1] = xd;
        for (i = 0; yn[N] > y[i]; i++);
        a[0] = y[i - 1], a[1] = y[i], a[2] = y[i + 1];
        b[0] = x[i - 1], b[1] = x[i], b[2] = x[i + 1];
        xn[N] = Lagrange(a, b, yn[N]);
        yn[N + 1] = R * xn[N] / (R + 1) + xd / (R + 1);
    }
    temp = yn[N];
    printf("精馏段理论板数：%d\n", N - 1);
    for (i = 1; i < 100; i++)
    {
        yn[i] = 0;
        xn[i] = 0;
    }
    //数组归零
    xq = dEt2(-xd / (R + 1), -1, xf / (q - 1), -1) / dEt2(R / (R + 1), -1, q / (q - 1), -1);
    yq = dEt2(R / (R + 1), -xd / (R + 1), q / (q - 1), xf / (q - 1)) / dEt2(R / (R + 1), -1, q / (q - 1), -1);
    for (N = 1; xn[N - 1] > xw; N++)    //逐板计算法求提馏段理论板数
    {
        yn[1] = temp;
        for (i = 0; yn[N] > y[i]; i++);
        a[0] = y[i - 1], a[1] = y[i], a[2] = y[i + 1];
        b[0] = x[i - 1], b[1] = x[i], b[2] = x[i + 1];
        xn[N] = Lagrange(a, b, yn[N]);
        yn[N + 1] = (xw - yq) / (xw - xq) * xn[N] - xq * (xw - yq) / (xw - xq) + yq;
    }
    printf("提馏段理论板数：%d\n", N - 1);
}
double Lagrange(double a[], double b[], double x)
{
    int i, j;
    double L = 0, l = 1;
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)//计算拉格朗日插值法每一项的基函数
        {
            if (i != j)
            {
                l *= (x - a[j]) / (a[i] - a[j]);
            }
        }
        L += b[i] * l;//每一项拉格朗日插值
        l = 1;
    }
    return L;
}
double dEt2(double a1, double b1,
            double a2, double b2)
{
    double d;
    d = a1 * b2 - b1 * a2;
    return d;
}
