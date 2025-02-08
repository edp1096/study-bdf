package main

import (
	"fmt"
	"math"
)

// 함수 f(x) 정의
func f(x float64) float64 {
	return math.Pow(x, 2) - 2
}

// 도함수 f'(x) 정의
func f_prime(x float64) float64 {
	return 2 * x
}

// 뉴턴-랩슨법 함수
func newton_raphson(x0 float64, tol float64, max_iter int) float64 {
	for i := 0; i < max_iter; i++ {
		x1 := x0 - f(x0)/f_prime(x0)

		fmt.Println("Mid:", x1, x0, tol)

		if math.Abs(x1-x0) < tol {
			return x1
		}
		x0 = x1
	}

	return math.NaN() // 수렴하지 않으면 NaN 반환
}

/*
함수:   f(x) = x² - 2
도함수: f'(x) = 2x
*/
func main() {
	// x0 := 1.5
	x0 := 2.0
	tol := 1e-6
	result := newton_raphson(x0, tol, 100)
	fmt.Println(result) // 약 1.4142135623730951 출력
}
