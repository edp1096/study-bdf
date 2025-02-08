커피의 온도 변화


* Formula
```r
T = 커피온도

dT/dt = -k(T-Troom)
T(t) = Troom + (T0 - Troom)e^(-kt)
```

* Forward Euler
```r
T[n+1] = T[n] - k*dt*(T[n] - Troom)
```

* Backward Euler
```r
T[n+1] = (T[n] + k*dt*Troom)/(1 + k*dt)
```

* Trapezoidal
```r
T[n+1] = (T[n]*(1 - k*dt/2) + k*dt*Troom)/(1 + k*dt/2)
```

* Gear's method (BDF)
```r
1차 (Backward Euler Method)
y^(n+1) = y^n + Δt f(y^(n+1))

2차
y^(n+1) =
(4/3)y^n
- (1/3)y^(n-1)
+ (2/3)Δt*f(y^(n+1))

3차
y^(n+1) =
(18/11)y^n
- (9/11)y^(n-1)
+ (2/11)y^(n-2)
+ (6/11)Δt*f(y^(n+1))

4차
y^(n+1) =
(48/25)y^n
- (36/25)y^(n-1)
+ (16/25)y^(n-2)
- (3/25)y^(n-3)
+ (12/25)Δt*f(y^(n+1))

5차
y^(n+1) =
(300/137)y^n
- (300/137)y^(n-1)
+ (200/137)y^(n-2)
- (75/137)y^(n-3)
+ (12/137)y^(n-4)
+ (60/137)Δt*f(y^(n+1))

6차
y^(n+1) =
(360/147)y^n - (450/147)y^(n-1)
+ (400/147)y^(n-2)
- (225/147)y^(n-3)
+ (72/147)y^(n-4)
- (10/147)y^(n-5)
+ (60/147)Δt*f(y^(n+1))


Extraolation Method:

1차
y[n+1] = y[n]

2차
y[n+1] = y[n] + (y[n] - y[n-1])

3차
y[n+1] = y[n] + (y[n] - y[n-1]) + (y[n] - 2y[n-1] + y[n-2])/2!

4차
y[n+1] = y[n] + (y[n] - y[n-1]) + (y[n] - 2y[n-1] + y[n-2])/2!
        + (y[n] - 3y[n-1] + 3y[n-2] - y[n-3])/3!

5차
y[n+1] = y[n] + (y[n] - y[n-1]) + (y[n] - 2y[n-1] + y[n-2])/2!
        + (y[n] - 3y[n-1] + 3y[n-2] - y[n-3])/3!
        + (y[n] - 4y[n-1] + 6y[n-2] - 4y[n-3] + y[n-4])/4!

6차
y[n+1] = y[n] + (y[n] - y[n-1]) + (y[n] - 2y[n-1] + y[n-2])/2!
        + (y[n] - 3y[n-1] + 3y[n-2] - y[n-3])/3!
        + (y[n] - 4y[n-1] + 6y[n-2] - 4y[n-3] + y[n-4])/4!
        + (y[n] - 5y[n-1] + 10y[n-2] - 10y[n-3] + 5y[n-4] - y[n-5])/5!
```
