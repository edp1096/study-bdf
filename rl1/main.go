package main

import (
	"fmt"
	"math"
)

// 회로 파라미터
const (
	R     = 100.0  // 저항 (Ω)
	L     = 1e-3   // 인덕턴스 (H)
	dt    = 0.0001 // 시간 간격 (s)
	tstop = 0.002  // 시뮬레이션 종료 시간 (s)
)

// 입력 전압 함수 (사인파)
func V_in(t float64) float64 {
	return 5 * math.Sin(2*math.Pi*1000*t) // 1kHz, 5V 입력 전압
}

// 전류 계산
func calculateCurrent(prevI, Vin, t float64) float64 {
	// 등가 전류원 방법
	// I(n+1) = I(n) + (Vin - R*I(n)) * dt/L
	return prevI + (Vin-R*prevI)*dt/L
}

// 뉴턴-랩슨
func newtonRaphsonCurrent(prevI, Vin, t float64) float64 {
	maxIter := 50
	tolerance := 1e-10

	I := prevI // 초기 추정값
	for iter := 0; iter < maxIter; iter++ {
		// 비선형 방정식의 함수값
		F := I - prevI - (Vin-R*I)*dt/L

		// 자코비안 (도함수)
		dF := 1.0 + R*dt/L

		// 뉴턴-랩슨 스텝
		delta := -F / dF

		// 새로운 추정값
		Inew := I + delta

		// 수렴 조건
		if math.Abs(delta) < tolerance {
			return Inew
		}

		I = Inew
	}

	return I
}

// 시뮬레이션 주요 함수
func simulateRLCircuit() ([]float64, []float64, []float64) {
	N := int(tstop / dt)
	t := make([]float64, N+1)
	I := make([]float64, N+1)
	Vin := make([]float64, N+1)
	VL := make([]float64, N+1)

	// 초기화
	I[0] = 0

	// 시간과 입력 전압 초기화
	for n := 0; n <= N; n++ {
		t[n] = float64(n) * dt
		Vin[n] = V_in(t[n])
	}

	// 각 시간 단계에서 전류 계산
	for n := 1; n <= N; n++ {
		// 뉴턴-랩슨 방법으로 전류 계산
		I[n] = newtonRaphsonCurrent(I[n-1], Vin[n], t[n])

		// 인덕터 전압 계산 (L * dI/dt)
		VL[n] = L * (I[n] - I[n-1]) / dt
	}

	return t, I, VL
}

func main() {
	fmt.Println("RL 회로 수치해석")
	fmt.Printf("저항: %.2f Ω\n", R)
	fmt.Printf("인덕턴스: %.3f H\n", L)
	fmt.Printf("시간 간격: %.4f s\n", dt)
	fmt.Printf("총 시뮬레이션 시간: %.4f s\n\n", tstop)

	t, I, VL := simulateRLCircuit()

	// 출력 (모든 데이터 포인트)
	fmt.Println(" 시간(s)   입력전압(V)   전류(A)   인덕터전압(V)")
	fmt.Println("----------------------------------------------")
	for n := 0; n < len(t); n += max(len(t)/20, 1) {
		fmt.Printf("%8.4f  %8.4f  %10.6f  %10.6f\n", t[n], V_in(t[n]), I[n], VL[n])
	}
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}
