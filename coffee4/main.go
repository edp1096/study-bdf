/*
사다리꼴법:
T[n+1] = (T[n]*(2 + k*dt) + 2*k*dt*Troom)/(2 + k*dt)
*/

package main

import (
	"fmt"
	"math"
)

func solveTrapezoidal(T0, Troom, k float64, dt float64, totalTime float64) []float64 {
	steps := int(totalTime/dt) + 1
	T := make([]float64, steps)
	T[0] = T0

	// Trapezoidal method
	// T[n+1] = (T[n]*(2 + k*dt) + 2*k*dt*Troom)/(2 + k*dt)
	for n := 0; n < steps-1; n++ {
		T[n+1] = (T[n]*(1-k*dt/2) + k*dt*Troom) / (1 + k*dt/2)
	}

	return T
}

func main() {
	T0 := 90.0        // 초기 온도
	Troom := 20.0     // 실내 온도
	k := 0.1          // 냉각 상수
	dt := 1.0         // 시간 간격
	totalTime := 30.0 // 총 시간

	// 해석해 계산
	analyticT := Troom + (T0-Troom)*math.Exp(-k*totalTime)

	fmt.Printf("총 단계 수: %d\n", int(totalTime/dt)+1)
	fmt.Printf("시간 간격: %.1f\n", dt)
	fmt.Printf("초기 온도: %.1f°C\n", T0)
	fmt.Printf("실내 온도: %.1f°C\n", Troom)
	fmt.Printf("냉각 상수: %.1f\n", k)
	fmt.Printf("총 시간: %.1f분\n\n", totalTime)

	fmt.Printf("시간(분)\t온도(°C)\n")
	solution := solveTrapezoidal(T0, Troom, k, dt, totalTime)
	for t := 0.0; t <= totalTime; t += dt {
		idx := int(t/dt + 0.5)
		fmt.Printf("%.1f\t%.2f\n", t, solution[idx])
	}

	finalT := solution[len(solution)-1]
	error := math.Abs(finalT - analyticT)
	relError := error / math.Abs(analyticT) * 100.0

	fmt.Printf("\n사다리꼴법 최종 결과: %.6f°C (err: %.6f%%)\n", finalT, relError)
	fmt.Printf("이론값: %.6f°C\n", analyticT)
}
