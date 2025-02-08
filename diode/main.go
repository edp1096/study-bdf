package main

import (
	"fmt"
	"math"
)

const (
	Is = 1e-12 // 역 포화 전류 (예시)
	n  = 1.5   // 이상계수 (예시)
	Vt = 0.026 // 열 전압 (실온에서)
)

func diodeCurrent(Vd float64) float64 {
	return Is * (math.Exp(Vd/(n*Vt)) - 1)
}

func diodeDynamicResistance(Vd float64) float64 {
	return Is / (n * Vt) * math.Exp(Vd/(n*Vt))
}

func main() {
	// Vd := 0.6 // 다이오드 전압
	Vd := 0.7 // 다이오드 전압
	current := diodeCurrent(Vd)
	resistance := diodeDynamicResistance(Vd)
	fmt.Printf("전류: %f A, 동적 저항: %f Ohm\n", current, resistance)
}
