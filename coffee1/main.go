package main

import (
	"fmt"
	"math"
)

func main() {
	T0 := 90.0    // 초기 온도
	Troom := 20.0 // 실내 온도
	k := 0.1      // 냉각 상수
	t := 30.0     // 시간(분)

	T := Troom + (T0-Troom)*math.Exp(-k*t) // 30분 후 온도 계산

	fmt.Printf("30분 후 온도: %.2f°C\n", T)
}
