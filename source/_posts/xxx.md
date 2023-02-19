---
title: xxx
date: 2023-02-19 23:09:17
tags:初等数论
---

# 初等数论 #

## 素数筛 ##

~~~c++
void sieve(int n)
{
    int m = (int)sqrt(n + 0.5);
    memset(vis, 0, sizeof(vis));
    for (int i = 2; i <= m; i++)
        if (!vis[i])
            for (int j = i * i;j <= n;j += i)
                vis[j] = 1;
}
int gen_primes(int n)
{
    sieve(n);
    int c = 0;
    for (int i = 2;i <= n;++i)
        if (!vis[i])
            prime[c++] = i;
    return c;
}//返回小于等于n的素数个数
~~~

## 扩展欧几里得算法 #

~~~c++
typedef long long ll;
void egcd(ll a, ll b, ll& d, ll& x, ll& y)
{
    if (!b) { d = a;x = 1;y = 0; }
    else { egcd(b, a % b, d, y, x);y -= x * (a / b); }
}//扩展gcd算法，d是最大公约数，ax+by=d，其中|x|+|y|最小
~~~

存在方程$ax+by=c$，且 $x,y,c$ 为整数
上述方程有解的充要条件是 $gcd(a,b)|c$ （裴蜀定理）,利用扩展欧几里得算法可以求得方程的解
若其中一组解为  $(x_0,y_0)$ ,它的任意整数解都可以写成 $(x_0+kb',y_0-ka')$ ,其中 $a'=a/gcd(a,b),b'=b/gcd(a,b)$

## 乘法逆元 ##

~~~c++
ll inv(ll a, ll n)
{
    ll d, x, y;
    egcd(a, n, d, x, y);
    return d == 1 ? (x + n) % n : -1;
}//计算模n下a的逆
~~~

若存在 $a^{-1}$ ，使得 $a*a^{-1}\equiv 1\pmod n$，则称 $a^{-1}$为 $a$ 模 $n$ 下的逆

## 欧拉函数 ##

$phi(i)= n(1-\frac 1{p_1})(1-\frac 1{p_2})...(1-\frac 1{p_k})$
其中， $p$ 为 $i$ 的质因数

~~~c++
int phi(int n)
{
    int m = (int)sqrt(n + 0.5);
    int ans = n;
    for (int i = 2;i <= m;++i)
        if (n % i == 0)
        {
            ans = ans / i * (i - 1);
            while (n % i == 0)
                n /= i;
        }
    if (n > 1)
        ans = ans / n * (n - 1);
}//欧拉函数等于不超过n且和n互素的整数个数
~~~

## 欧拉定理求乘法逆元 ##

~~~c++
typedef long long ll;
ll eulerinv(ll a, ll n)
{
    ll t = pow(a, phi(n) - 1);
    return t % n;
}//a与n互质，计算模n下a的逆
~~~

欧拉定理：若a与n互质，则 $a^{phi(n)}\equiv 1\pmod n$

## 用类似素数筛的方法计算欧拉函数表 ##

~~~c++
void phi_table(int n)
{
    phi[1] = 1;
    for (int i = 2;i <= n;++i)
        if (!phi[i])
            for (int j = i;j <= n;j += i)
            {
                if (!phi[j])
                    phi[j] = j;
                phi[j] = phi[j] / i * (i - 1);
            }
}
~~~

## 算术基本定理 ##

1).对于任意的大于1的自然数N，都有 N = p1^a1 * p2^a2 * p3^a3……pn^an(p1,……
pn代表素数,a1,a2……an代表N中ai的个数)，列如：12=2^2 *3^1,p1=2,a1=2,p2=3
,a2=1。
2）N的正约数个数=(1+a1)*(1+a2)*(1+a3)*……*（1+an）.
3）N的所有正约数和为:(1+p1+p1^2+……+p1^c1)*(1+p2+p2^2+……+p2^c2)*……*
(1+pn+pn^2+……+pn^cn).