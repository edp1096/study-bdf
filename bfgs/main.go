package main

import (
	"fmt"
	"math"
)

// f(x) = x^2 - 2
func f(x float64) float64 {
	return x*x - 2
}

// f'(x) = 2x
func df(x float64) float64 {
	return 2 * x
}

// 벡터 내적
func dot(a, b []float64) float64 {
	result := 0.0
	for i := range a {
		result += a[i] * b[i]
	}
	return result
}

// 벡터의 크기(norm)
func norm(v []float64) float64 {
	return math.Sqrt(dot(v, v))
}

// 벡터 덧셈
func add(a, b []float64) []float64 {
	result := make([]float64, len(a))
	for i := range a {
		result[i] = a[i] + b[i]
	}
	return result
}

// 벡터 뺄셈
func subtract(a, b []float64) []float64 {
	result := make([]float64, len(a))
	for i := range a {
		result[i] = a[i] - b[i]
	}
	return result
}

// 벡터 스칼라 곱
func scale(scalar float64, v []float64) []float64 {
	result := make([]float64, len(v))
	for i := range v {
		result[i] = scalar * v[i]
	}
	return result
}

// 단순화된 line search
func lineSearch(x []float64, direction []float64) float64 {
	alpha := 1.0
	beta := 0.5
	c := 1e-4

	fx := f(x[0])
	gradient := []float64{df(x[0])}
	init_slope := dot(gradient, direction)

	for i := 0; i < 20; i++ { // 최대 20번 시도
		new_x := add(x, scale(alpha, direction))
		if f(new_x[0]) <= fx+c*alpha*init_slope {
			return alpha
		}
		alpha *= beta
	}
	return alpha
}

// BFGS 업데이트
func updateHessian(B [][]float64, s, y []float64) [][]float64 {
	n := len(B)
	newB := make([][]float64, n)
	for i := range newB {
		newB[i] = make([]float64, n)
		for j := range newB[i] {
			newB[i][j] = B[i][j]
		}
	}

	sy := dot(s, y)
	if sy <= 0 {
		return B // 업데이트 건너뛰기
	}

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			term1 := (sy + dot(y, y)) * s[i] * s[j] / (sy * sy)
			term2 := (B[i][0]*s[0]*y[j] + y[i]*s[0]*B[0][j]) / sy
			newB[i][j] = B[i][j] + term1 - term2
		}
	}
	return newB
}

// BFGS 알고리즘
func bfgs(x0 float64, tol float64, maxIter int) float64 {
	x := []float64{x0}
	// n := 1

	// 초기 Hessian 근사
	B := [][]float64{{1.0}}

	for iter := 0; iter < maxIter; iter++ {
		grad := []float64{df(x[0])}
		if norm(grad) < tol {
			fmt.Printf("Converged after %d iterations\n", iter+1)
			return x[0]
		}

		direction := scale(-1.0/B[0][0], grad) // 탐색 방향 계산
		alpha := lineSearch(x, direction)      // line search로 step size 결정

		// 새로운 위치로 이동
		xNew := add(x, scale(alpha, direction))
		gradNew := []float64{df(xNew[0])}

		// s와 y 벡터 계산
		s := subtract(xNew, x)
		y := subtract(gradNew, grad)

		B = updateHessian(B, s, y) // Hessian 업데이트
		x = xNew                   // 현재 위치 업데이트

		fmt.Printf("Iteration %d: x = %v, f(x) = %v\n", iter+1, x[0], f(x[0]))
	}

	fmt.Println("Maximum iterations reached")
	return x[0]
}

func main() {
	x0 := 2.0      // 초기값
	tol := 1e-6    // 허용오차
	maxIter := 100 // 최대 반복 횟수

	result := bfgs(x0, tol, maxIter)
	fmt.Printf("Final result: x = %v\n", result)
	fmt.Printf("f(x) = %v\n", f(result))
	fmt.Printf("Expected sqrt(2) = %v\n", math.Sqrt(2))
}
