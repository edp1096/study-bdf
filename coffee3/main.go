/*
후방오일러:
T[n+1] = (T[n] + k*dt*Troom)/(1 + k*dt)
*/

package main

import (
	"fmt"
)

func main() {
	T := 90.0     // 초기 온도 (°C)
	Troom := 20.0 // 실내 온도 (°C)
	k := 0.1      // 냉각 속도 상수

	dt := 1.0    // 시간 간격 (분)
	tEnd := 30.0 // 총 시뮬레이션 시간 (분)

	t := 0.0
	fmt.Printf("시간(분)\t온도(°C)\n")
	fmt.Printf("%.1f\t%.2f\n", t, T)

	for t < tEnd {
		// Backward Euler method
		T = (T + k*dt*Troom) / (1 + k*dt)
		t = t + dt

		fmt.Printf("%.1f\t%.2f\n", t, T)
	}
}
