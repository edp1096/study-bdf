/*
https://en.wikipedia.org/wiki/Backward_differentiation_formula

Gear's Method: BDF + Newton-Raphson Method

BDF(Backward Differentiation Formula):

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
*/

package main

import (
	"fmt"
	"math"
)

type bdfMethod struct {
	coefficients []float64 // BDF 계수
	beta         float64   // BDF의 f(y^(n+1)) 계수
}

// BDF 차수별 계수
var bdfMethods = []bdfMethod{
	{[]float64{1.0}, 1.0},                                                                                                // 1차
	{[]float64{4.0 / 3.0, -1.0 / 3.0}, 2.0 / 3.0},                                                                        // 2차
	{[]float64{18.0 / 11.0, -9.0 / 11.0, 2.0 / 11.0}, 6.0 / 11.0},                                                        // 3차
	{[]float64{48.0 / 25.0, -36.0 / 25.0, 16.0 / 25.0, -3.0 / 25.0}, 12.0 / 25.0},                                        // 4차
	{[]float64{300.0 / 137.0, -300.0 / 137.0, 200.0 / 137.0, -75.0 / 137.0, 12.0 / 137.0}, 60.0 / 137.0},                 // 5차
	{[]float64{360.0 / 147.0, -450.0 / 147.0, 400.0 / 147.0, -225.0 / 147.0, 72.0 / 147.0, -10.0 / 147.0}, 60.0 / 147.0}, // 6차
}

// 미분방정식: dT/dt = -k(T - Troom)
func df(_ float64, T, k, Troom float64) float64 {
	return -k * (T - Troom)
}

// 4th order runge-kutta
func rk4Step(t, T, dt, k, Troom float64) float64 {
	k1 := df(t, T, k, Troom)
	k2 := df(t+dt/2, T+dt*k1/2, k, Troom)
	k3 := df(t+dt/2, T+dt*k2/2, k, Troom)
	k4 := df(t+dt, T+dt*k3, k, Troom)

	return T + dt*(k1+2*k2+2*k3+k4)/6
}

// 외삽법 - Newton's divided differences
func extrapolate(T []float64, n, order int) float64 {
	if n < order-1 {
		return T[n] // 초기에는 이전 값 사용
	}

	// 1차 차분
	diffs := make([][]float64, order)
	diffs[0] = make([]float64, order)
	for i := 0; i < order; i++ {
		diffs[0][i] = T[n-i]
	}

	// 고차 차분 계산
	for i := 1; i < order; i++ {
		diffs[i] = make([]float64, order-i)
		for j := 0; j < order-i; j++ {
			diffs[i][j] = diffs[i-1][j] - diffs[i-1][j+1]
		}
	}

	// 외삽값
	factorial := 1.0
	pred := T[n]
	for i := 1; i < order; i++ {
		factorial *= float64(i)
		pred += diffs[i][0] / factorial
	}

	return pred
}

func solveBDF(order int, T0, Troom, k float64, dt float64, totalTime float64) []float64 {
	steps := int(totalTime/dt) + 1
	T := make([]float64, steps)
	T[0] = T0

	// RK4로 초기값 준비
	for i := 1; i < order; i++ {
		T[i] = rk4Step(float64(i-1)*dt, T[i-1], dt, k, Troom)
	}

	method := bdfMethods[order-1]

	// BDF
	for n := order - 1; n < steps-1; n++ {
		// Predictor - 뉴턴분할차분법
		Tn := extrapolate(T, n, order)

		// Corrector - BDF + 뉴턴-랩슨
		maxIter := 50
		tolerance := 1e-10

		for iter := 0; iter < maxIter; iter++ {
			// 이전 값들의 선형 조합 계산
			sum := 0.0
			for i := 0; i < len(method.coefficients); i++ {
				sum += method.coefficients[i] * T[n-i]
			}

			// 현재 시점의 도함수 값
			fn := df(float64(n+1)*dt, Tn, k, Troom)
			rightSide := sum + dt*method.beta*fn

			F := Tn - rightSide          // 잔차(residual)
			dF := 1.0 + dt*method.beta*k // 자코비안
			delta := F / dF              // 뉴턴 스텝
			Tnew := Tn - delta           // 업데이트

			// 수렴 확인
			if math.Abs(Tnew-Tn) < tolerance {
				Tn = Tnew
				break
			}

			// 발산 체크
			if math.IsNaN(Tnew) || math.Abs(Tnew-Troom) > 10*math.Abs(T0-Troom) {
				Tn = T[n]
				break
			}

			Tn = Tnew
		}
		T[n+1] = Tn
	}
	return T
}

func main() {
	T0 := 90.0        // 초기 온도
	Troom := 20.0     // 실내 온도
	k := 0.1          // 냉각 상수
	dt := 1.0         // 시간 간격
	totalTime := 30.0 // 총 시간

	analyticT := Troom + (T0-Troom)*math.Exp(-k*totalTime)

	fmt.Printf("총 단계 수: %d\n", int(totalTime/dt)+1)
	fmt.Printf("시간 간격: %.3f\n", dt)
	fmt.Printf("초기 온도: %.1f°C\n", T0)
	fmt.Printf("실내 온도: %.1f°C\n", Troom)
	fmt.Printf("냉각 상수: %.1f\n", k)
	fmt.Printf("총 시간: %.1f분\n\n", totalTime)

	for order := 1; order <= 6; order++ {
		solution := solveBDF(order, T0, Troom, k, dt, totalTime)
		finalT := solution[len(solution)-1]

		error := math.Abs(finalT - analyticT)
		relError := error / math.Abs(analyticT) * 100.0
		fmt.Printf("Gear%d: %.6f°C (err: %.6f%%)\n", order, finalT, relError)
	}

	fmt.Println()
	fmt.Printf("이론값: %.6f°C\n", analyticT)
}
