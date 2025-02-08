/*
전방오일러:
T[n+1] = T[n] - k*dt*(T[n] - Troom)
*/

package main

import (
	"fmt"
	"math"
)

func dTdt(T, Troom, k float64) float64 {
	return -k * (T - Troom)
}

func main() {
	// 초기 조건
	T := 90.0     // 초기 온도 (°C)
	Troom := 20.0 // 실내 온도 (°C)
	k := 0.1      // 냉각 속도 상수

	dt := 0.1    // 시간 간격 (분)
	tEnd := 30.0 // 총 시뮬레이션 시간 (분)

	t := 0.0
	fmt.Printf("시간(분)\t온도(°C)\t변화율(°C/분)\n")
	fmt.Printf("%.1f\t%.2f\t%.2f\n", t, T, dTdt(T, Troom, k))

	for t < tEnd {
		// Forward Euler method
		dT := dTdt(T, Troom, k)
		T = T + dt*dT
		t = t + dt

		// 1분 단위로 결과 출력
		if math.Mod(t, 1.0) < dt/2 {
			fmt.Printf("%.1f\t%.2f\t%.2f\n", t, T, dT)
		}
	}
}
