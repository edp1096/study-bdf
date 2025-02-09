package main

import (
	"fmt"
	"math"
)

// 회로 파라미터
const (
	R = 100.0 // 저항 (Ω)
	L = 1e-3  // 인덕턴스 (H)
	// dt    = 1e-6   // 타임스텝 (s)
	dt    = 1e-5   // 타임스텝 (s)
	tstop = 0.002  // 종료 시간 (s)
	Vpeak = 5.0    // 전압 피크 (V)
	freq  = 1000.0 // 주파수 (Hz)
)

// BDF 계수 구조체
type bdfMethod struct {
	coefficients []float64
	beta         float64
}

var bdfMethods = []bdfMethod{
	{[]float64{1.0}, 1.0},                                                                                                // 1차
	{[]float64{4.0 / 3.0, -1.0 / 3.0}, 2.0 / 3.0},                                                                        // 2차
	{[]float64{18.0 / 11.0, -9.0 / 11.0, 2.0 / 11.0}, 6.0 / 11.0},                                                        // 3차
	{[]float64{48.0 / 25.0, -36.0 / 25.0, 16.0 / 25.0, -3.0 / 25.0}, 12.0 / 25.0},                                        // 4차
	{[]float64{300.0 / 137.0, -300.0 / 137.0, 200.0 / 137.0, -75.0 / 137.0, 12.0 / 137.0}, 60.0 / 137.0},                 // 5차
	{[]float64{360.0 / 147.0, -450.0 / 147.0, 400.0 / 147.0, -225.0 / 147.0, 72.0 / 147.0, -10.0 / 147.0}, 60.0 / 147.0}, // 6차
}

// 전압원
func Vin(t float64) float64 {
	return Vpeak * math.Sin(2.0*math.Pi*freq*t)
}

// 저항 전류 계산
func calcResistorCurrent(v float64) float64 {
	return v / R
}

// 인덕터 전압 계산
func calcInductorVoltage(di, dt float64) float64 {
	return L * di / dt
}

// Trapezoidal 적분을 사용한 초기 전류 계산
func calcInitialCurrent(v0, v1, i0 float64, dt float64) float64 {
	// v = L*di/dt + Ri 방정식에서
	// i1 = i0 + (dt/2L)*(v0 + v1 - R*(i0 + i1))
	// => i1 = (i0 + (dt/2L)*(v0 + v1 - R*i0))/(1 + R*dt/2L)
	return (i0 + (dt/(2*L))*(v0+v1-R*i0)) / (1 + R*dt/(2*L))
}

// 상태 벡터 구조체
type State struct {
	t     float64 // 시간
	v     float64 // 전압
	i     float64 // 전류
	di_dt float64 // 전류 변화율
}

func newState(t float64) State {
	return State{t: t}
}

func solveBDF(order int) []State {
	steps := int(tstop/dt) + 1
	states := make([]State, steps)

	// DC 동작점 계산
	states[0] = newState(0)
	states[0].v = Vin(0)
	states[0].i = calcResistorCurrent(states[0].v)
	states[0].di_dt = 0

	// Trapezoidal로 첫 2스텝 계산
	for n := 0; n < 2; n++ {
		states[n+1] = newState(float64(n+1) * dt)
		states[n+1].v = Vin(states[n+1].t)
		states[n+1].i = calcInitialCurrent(
			states[n].v,
			states[n+1].v,
			states[n].i,
			dt,
		)
		states[n+1].di_dt = (states[n+1].i - states[n].i) / dt
	}

	// BDF 메인 루프
	history := make([]float64, 6) // 이전 전류값 저장용

	for n := 2; n < steps-1; n++ {
		states[n+1] = newState(float64(n+1) * dt)
		states[n+1].v = Vin(states[n+1].t)

		// 현재 BDF 차수 결정
		currentOrder := order
		if n < order {
			currentOrder = n
		}
		method := bdfMethods[currentOrder-1]

		// 이전 전류값 수집
		for i := 0; i < currentOrder; i++ {
			history[i] = states[n-i].i
		}

		// BDF로 새로운 전류값 계산
		// Σ(α[k]*i[n-k]) + β*dt/L*(v[n+1] - R*i[n+1]) = i[n+1]
		sum := 0.0
		for i := 0; i < len(method.coefficients); i++ {
			sum += method.coefficients[i] * history[i]
		}

		// i[n+1] = (sum + β*dt*v[n+1]/L)/(1 + β*dt*R/L)
		states[n+1].i = (sum + method.beta*dt*states[n+1].v/L) /
			(1 + method.beta*dt*R/L)

		// 전류 변화율과 인덕터 전압 계산
		states[n+1].di_dt = (states[n+1].i - states[n].i) / dt

		// LTE 계산 및 안정성 검사는 여기서 추가 가능
	}

	return states
}

func main() {
	fmt.Println("RL 회로 수치해석 (BDF/Gear Method)")
	fmt.Printf("R = %.1f Ω\n", R)
	fmt.Printf("L = %.3f H\n", L)
	fmt.Printf("Timestep = %.6f s\n", dt)
	fmt.Printf("End time = %.3f s\n\n", tstop)

	// 이론적 최대값 계산
	maxDi_dt := Vpeak * 2 * math.Pi * freq / math.Sqrt(R*R+(2*math.Pi*freq*L)*(2*math.Pi*freq*L))
	theoreticalMaxVL := L * maxDi_dt
	fmt.Printf("이론적 최대 인덕터 전압 = %.6f V\n\n", theoreticalMaxVL)

	for order := 1; order <= 6; order++ {
		fmt.Printf("Gear%d:\n", order)

		states := solveBDF(order)

		// 결과 분석
		maxVL := 0.0
		for i := 1; i < len(states); i++ {
			vL := L * states[i].di_dt
			if math.Abs(vL) > math.Abs(maxVL) {
				maxVL = vL
			}
		}

		fmt.Printf("최대 인덕터 전압: %.6f V\n", maxVL)
		fmt.Printf("오차율: %.6f%%\n\n", (math.Abs(maxVL)-theoreticalMaxVL)/theoreticalMaxVL*100)

		// 결과 출력
		fmt.Println("   시간(s)    전압(V)    전류(A)    di/dt(A/s)    vL(V)")
		fmt.Println("---------------------------------------------------------")

		printStep := int(math.Max(float64(len(states)/20), 1))
		for i := 0; i < len(states); i += printStep {
			s := states[i]
			vL := L * s.di_dt
			fmt.Printf("%10.6f %10.6f %10.6f %10.2f %10.6f\n",
				s.t, s.v, s.i, s.di_dt, vL)
		}
		fmt.Println()
	}
}
