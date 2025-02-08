package main

import (
	"fmt"
	"math"
)

// 회로 파라미터
const (
	R     = 100.0   // 저항 (Ω)
	L     = 1e-3    // 인덕턴스 (H)
	dt    = 0.00002 // 타임스텝 (s)
	tstop = 0.002   // 시뮬레이션 종료 시간 (s)
)

type bdfMethod struct {
	coefficients []float64 // BDF 계수
	beta         float64   // f(y^(n+1)) 계수
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

// 입력 전압 함수 (사인파)
func V_in(t float64) float64 {
	return 5 * math.Sin(2*math.Pi*1000*t) // 1kHz, 5V 입력 전압
}

// RL 회로 미분방정식
func dI(t, I float64) float64 {
	return (V_in(t) - R*I) / L
}

// RK4 초기값 준비 메서드 (개선된 버전)
func rk4Step(t, I, dt float64) float64 {
	k1 := dI(t, I)
	k2 := dI(t+dt/2, I+dt*k1/2)
	k3 := dI(t+dt/2, I+dt*k2/2)
	k4 := dI(t+dt, I+dt*k3)

	return I + dt*(k1+2*k2+2*k3+k4)/6
}

// 개선된 외삽법 (안정성 강화)
func extrapolate(values []float64, n, order int) float64 {
	if n < order-1 {
		return values[n]
	}

	// 가중치 기반 외삽
	pred := values[n]
	weight := 1.0
	factorial := 1.0

	for i := 1; i < order; i++ {
		factorial *= float64(i)
		diff := values[n] - values[n-i]
		pred += diff / (factorial * math.Pow(float64(i), 1.5))
		weight += 1.0 / (factorial * math.Pow(float64(i), 1.5))
	}

	return pred / weight
}

// BDF 방법으로 RL 회로 해석 (개선된 버전)
func solveBDF(order int) ([]float64, []float64, []float64, []float64) {
	N := int(math.Floor(tstop/dt)) + 1
	t := make([]float64, N+1)
	I := make([]float64, N+1)
	V_L := make([]float64, N+1)
	V_in_arr := make([]float64, N+1)

	// 시간 배열 초기화
	for n := 0; n <= N; n++ {
		t[n] = float64(n) * dt
		V_in_arr[n] = V_in(t[n])
	}

	// 초기값 설정 개선
	I[0] = 0 // 명시적 초기 전류 0 설정

	// RK4로 초기값 준비 (각 차수에 맞게 조정)
	for i := 1; i < order; i++ {
		I[i] = rk4Step(t[i-1], I[i-1], dt)
	}

	method := bdfMethods[order-1]

	// BDF 방법 적용 (수렴성 및 안정성 강화)
	for n := order - 1; n < N; n++ {
		// 초기 추정값 (개선된 외삽법)
		I_pred := extrapolate(I, n, order)

		maxIter := 200
		tolerance := 1e-14
		dampingFactor := 1.0

		for iter := 0; iter < maxIter; iter++ {
			// 이전 값들의 선형 조합
			sum := 0.0
			for i := 0; i < len(method.coefficients); i++ {
				sum += method.coefficients[i] * I[n-i]
			}

			// 현재 시점 미분값
			fn := dI(t[n+1], I_pred)
			rightSide := sum + dt*method.beta*fn

			F := I_pred - rightSide                                     // 개선된 잔차 계산
			dF := 1.0 + dt*method.beta*R/L*(1.0+math.Abs(I_pred)/100.0) // 동적 자코비안 근사
			delta := -F / dF * dampingFactor                            // 뉴턴-랩슨 스텝
			I_new := I_pred + delta                                     // 전류값 갱신

			// 수렴 및 안정성 확인
			if math.Abs(delta) < tolerance {
				I_pred = I_new
				break
			}

			// 발산 체크 및 방지
			if iter > maxIter/2 {
				dampingFactor *= 0.5
			}

			if math.IsNaN(I_new) || math.Abs(I_new) > 1000 {
				I_pred = I[n]
				break
			}

			I_pred = I_new
		}

		I[n+1] = I_pred

		// 인덕터 전압 계산
		V_L[n+1] = L * (I[n+1] - I[n]) / dt
	}

	return t, I, V_L, V_in_arr
}

func main() {
	fmt.Println("RL 회로 BDF 수치해석")
	fmt.Printf("저항: %.2f Ω\n", R)
	fmt.Printf("인덕턴스: %.3f H\n", L)
	fmt.Printf("시간 간격: %.4f s\n", dt)
	fmt.Printf("총 시뮬레이션 시간: %.4f s\n\n", tstop)

	// 이론적 최대 인덕터 전압 계산
	maxDI := 2 * math.Pi * 1000 * 5 / R
	theoreticalMaxVL := L * maxDI

	// 각 BDF 차수별 결과 분석
	for order := 1; order <= 6; order++ {
		fmt.Printf("Gear%d:\n", order)
		t, I, V_L, V_in_arr := solveBDF(order)

		// 최대 인덕터 전압 계산
		maxVL := 0.0
		for _, vl := range V_L {
			if math.Abs(vl) > math.Abs(maxVL) {
				maxVL = vl
			}
		}

		// 오차율 계산
		error := math.Abs(maxVL) - math.Abs(theoreticalMaxVL)
		relError := error / math.Abs(theoreticalMaxVL) * 100.0

		fmt.Printf("최대 인덕터 전압: %.6f V\n", maxVL)
		fmt.Printf("최대 인덕터 전압 이론치: %.6f V\n", theoreticalMaxVL)
		fmt.Printf("오차율: %.6f%%\n\n", relError)

		printStep := max(len(t)/20, 1)
		fmt.Println(" 시간(s)   입력전압(V)   전류(A)   인덕터전압(V)")
		fmt.Println("----------------------------------------------")
		for n := 0; n < len(t); n += printStep {
			fmt.Printf("%8.4f  %8.4f  %10.6f  %10.6f\n", t[n], V_in_arr[n], I[n], V_L[n])
		}
		fmt.Println()
	}
}
