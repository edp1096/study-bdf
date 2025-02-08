package main

import (
	"fmt"
	"math"
)

func sqrt2TaylorSeries(terms int) float64 {
	result := 0.0
	factor := 1.0

	for n := 0; n < terms; n++ {
		if n > 0 {
			factor *= (0.5 - float64(n-1)) / float64(n) // Update the factorial
		}
		result += factor
	}

	return result
}

func main() {
	terms := 200 // Terms in the series count

	approxSqrt2 := sqrt2TaylorSeries(terms)

	fmt.Printf("Approximated sqrt(2) with %d terms: %.10f\n", terms, approxSqrt2)
	fmt.Printf("Actual sqrt(2): %.10f\n", math.Sqrt(2))
}
