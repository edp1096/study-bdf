package main

import (
	"fmt"
	"math"
)

func f(x float64) float64 {
	return math.Sin(x) - 0.5*x
}

func f_prime(x float64) float64 {
	return math.Cos(x) - 0.5
}

func newton_raphson(x0 float64, tol float64, max_iter int) float64 {
	for i := 0; i < max_iter; i++ {
		x1 := x0 - f(x0)/f_prime(x0)

		fmt.Println("Mid:", x1, x0, tol)

		if math.Abs(x1-x0) < tol {
			return x1
		}
		x0 = x1
	}
	return math.NaN()
}

/*
함수:   f(x) = sin(x) - x/2
도함수: f'(x) = cos(x) - 0.5
*/
func main() {
	x0 := 2.0 // 초기값 설정
	tol := 1e-6
	result := newton_raphson(x0, tol, 100)
	fmt.Println(result)
}
