package main

import (
	"fmt"
	"math"
)

const (
	R     = 100.0   // 저항 (Ω)
	L     = 1e-3    // 인덕턴스 (H)
	dt    = 0.00001 // 타임스텝 (s)
	tstop = 0.002   // 시뮬레이션 종료 시간 (s)
)

type bdfMethod struct {
	coefficients []float64
	beta         float64
}

var bdfMethods = []bdfMethod{
	{[]float64{1.0}, 1.0},                                                                                                // 1차
	{[]float64{4.0 / 3.0, -1.0 / 3.0}, 2.0 / 3.0},                                                                        // 2차
	{[]float64{18.0 / 11.0, -9.0 / 11.0, 2.0 / 11.0}, 6.0 / 11.0},                                                        // 3차
	{[]float64{48.0 / 25.0, -36.0 / 25.0, 16.0 / 25.0, -3.0 / 25.0}, 12.0 / 25.0},                                        // 4차
	{[]float64{300.0 / 137.0, -300.0 / 137.0, 200.0 / 137.0, -75.0 / 137.0, 12.0 / 137.0}, 60.0 / 137.0},                 // 5차
	{[]float64{360.0 / 147.0, -450.0 / 147.0, 400.0 / 147.0, -225.0 / 147.0, 72.0 / 147.0, -10.0 / 147.0}, 60.0 / 147.0}, // 6차
}

func V_in(t float64) float64 {
	return 5 * math.Sin(2*math.Pi*1000*t)
}

// 미분방정식
func dI(t, I float64) float64 {
	return (V_in(t) - R*I) / L
}

func rk4Step(t, I, dt float64) float64 {
	k1 := dI(t, I)
	k2 := dI(t+dt/2, I+dt*k1/2)
	k3 := dI(t+dt/2, I+dt*k2/2)
	k4 := dI(t+dt, I+dt*k3)
	return I + dt*(k1+2*k2+2*k3+k4)/6
}

// Newton's divided differences를 이용한 외삽법
func extrapolate(I []float64, n, order int) float64 {
	if n < order-1 {
		return I[n]
	}

	diffs := make([][]float64, order)
	diffs[0] = make([]float64, order)
	for i := 0; i < order; i++ {
		diffs[0][i] = I[n-i]
	}

	for i := 1; i < order; i++ {
		diffs[i] = make([]float64, order-i)
		for j := 0; j < order-i; j++ {
			diffs[i][j] = diffs[i-1][j] - diffs[i-1][j+1]
		}
	}

	factorial := 1.0
	pred := I[n]
	for i := 1; i < order; i++ {
		factorial *= float64(i)
		pred += diffs[i][0] / factorial
	}
	return pred
}

func solveBDF(targetOrder int) ([]float64, []float64, []float64) {
	steps := int(math.Floor(tstop/dt)) + 1
	t := make([]float64, steps)
	I := make([]float64, steps)
	V_L := make([]float64, steps)

	for n := 0; n < steps; n++ {
		t[n] = float64(n) * dt
	}

	// 초기값 - Trapezoidal
	I[1] = I[0] + (dt/(2*L))*(V_in(t[0])+V_in(t[1])-R*I[0])
	V_L[1] = L * (I[1] - I[0]) / dt

	// // 초기값 - RK4
	// for i := 1; i < targetOrder; i++ {
	// 	I[i] = rk4Step(t[i-1], I[i-1], dt)
	// 	V_L[i] = L * (I[i] - I[i-1]) / dt
	// }

	// BDF 메인 루프
	for n := 1; n < steps-1; n++ {
		currentOrder := targetOrder
		if n < targetOrder-1 {
			currentOrder = n + 1
		}

		method := bdfMethods[currentOrder-1]

		// Predictor - 외삽법으로 초기값 예측
		I_pred := extrapolate(I, n, currentOrder)

		// Corrector - BDF + 뉴턴-랩슨
		maxIter := 50
		tolerance := 1e-10
		// v_next := V_in(t[n+1])

		for iter := 0; iter < maxIter; iter++ {
			// 이전 값들의 선형 조합
			sum := 0.0
			for i := 0; i < len(method.coefficients); i++ {
				sum += method.coefficients[i] * I[n-i]
			}

			fn := dI(t[n+1], I_pred)
			rightSide := sum + dt*method.beta*fn
			F := I_pred - rightSide
			dF := 1.0 + dt*method.beta*R/L
			delta := -F / dF
			I_new := I_pred + delta

			// 수렴 확인
			if math.Abs(delta) < tolerance {
				I_pred = I_new
				break
			}

			// 발산 체크 및 방지
			maxAllowedCurrent := 5.0 / R * 10 // 이론적 최대전류의 10배
			if math.IsNaN(I_new) || math.Abs(I_new) > maxAllowedCurrent {
				I_pred = I[n] // 이전 값으로 복귀
				break
			}

			I_pred = I_new
		}

		I[n+1] = I_pred
		// 전압은 단순 차분으로 계산
		V_L[n+1] = L * (I[n+1] - I[n]) / dt
	}

	return t, I, V_L
}

func main() {
	fmt.Println("RL 회로 BDF 수치해석")
	fmt.Printf("저항: %.2f Ω\n", R)
	fmt.Printf("인덕턴스: %.3f H\n", L)
	fmt.Printf("시간 간격: %.4f s\n", dt)
	fmt.Printf("총 시뮬레이션 시간: %.4f s\n\n", tstop)

	maxDI := 2 * math.Pi * 1000 * 5 / R
	theoreticalMaxVL := L * maxDI

	for order := 1; order <= 6; order++ {
		fmt.Printf("Gear%d:\n", order)
		t, I, V_L := solveBDF(order)

		maxVL := 0.0
		for _, vl := range V_L {
			if math.Abs(vl) > math.Abs(maxVL) {
				maxVL = vl
			}
		}

		errRate := math.Abs(maxVL) - math.Abs(theoreticalMaxVL)
		relError := errRate / math.Abs(theoreticalMaxVL) * 100.0

		fmt.Printf("최대 인덕터 전압: %.6f V\n", maxVL)
		fmt.Printf("최대 인덕터 전압 이론치: %.6f V\n", theoreticalMaxVL)
		fmt.Printf("err: %.6f%%\n\n", relError)

		printStep := int(math.Max(float64(len(t)/20), 1))
		fmt.Println(" 시간(s)   입력전압(V)   전류(A)   V_L(V)")
		fmt.Println("----------------------------------------------")
		for n := 0; n < len(t); n += printStep {
			fmt.Printf("%8.4f  %8.4f  %10.6f  %7.3f\n", t[n], V_in(t[n]), I[n], V_L[n])
		}
		fmt.Println()
	}
}
