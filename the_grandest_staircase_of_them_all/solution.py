def solution(n):
	m = [[0 for i in range(n + 1)] for j in range(n + 1)]
	m[0][0] = 1  # base case
	
	for last in range(1, n + 1):
		for left in range(0, n + 1):
			m[last][left] = m[last - 1][left]
			if left >= last:
				m[last][left] += m[last - 1][left - last]
	          	
	return (m[n][n] -1)

print(solution(5))