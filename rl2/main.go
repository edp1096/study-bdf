package main

import (
	"fmt"
	"math"
)

// RL 회로 값 설정
const (
	R     = 100.0   // 저항 (Ohm)
	L     = 1e-3    // 인덕턴스 (H)
	dt    = 0.00001 // 타임스텝 (s)
	tstop = 0.002   // 시뮬레이션 종료 시간 (s)
)

func V_in(t float64) float64 {
	return 5 * math.Sin(2*math.Pi*1000*t) // 1kHz, 5V 입력 전압
}

// 자코비안 계산 함수
func computeJacobian(R, L, dt float64) float64 {
	return 1.0 + (2.0/3.0)*(dt/L)*R
}

func bdf2RLSimulation(R, L, dt, tstop float64) ([]float64, []float64, []float64) {
	N := int(tstop / dt)
	t := make([]float64, N+1)
	i := make([]float64, N+1)
	V_L := make([]float64, N+1)

	// 시간 배열 초기화
	for n := 0; n <= N; n++ {
		t[n] = float64(n) * dt
	}

	// 첫 번째 스텝 (Trapezoidal Method로 초기화)
	i[1] = i[0] + (dt/(2*L))*(V_in(t[0])+V_in(t[1])-R*i[0])

	// BDF 2차 반복
	for n := 1; n < N; n++ {
		v_next := V_in(t[n+1])
		i_next := i[n]

		for k := 0; k < 10; k++ { // Newton-Raphson 반복
			F := i_next - (4.0/3.0)*i[n] + (1.0/3.0)*i[n-1] - (2.0/3.0)*(dt/L)*(v_next-R*i_next)
			dF_di := computeJacobian(R, L, dt)
			delta := -F / dF_di
			i_next += delta
			if math.Abs(delta) < 1e-6 {
				break
			}
		}
		i[n+1] = i_next
		V_L[n+1] = L * (i[n+1] - i[n]) / dt
	}

	return t, i, V_L
}

func main() {
	t, i, V_L := bdf2RLSimulation(R, L, dt, tstop)

	for n := 0; n < len(t); n++ {
		fmt.Printf("time = %.4f s, current = %.6f A, V_L = %.6f V\n", t[n], i[n], V_L[n])
	}
}
