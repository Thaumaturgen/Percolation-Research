def binomialcoefficient (n,k):
	if k < 0 or k > n: 
		return 0
	if k > n -k :
		k = n-k 
	c=1
	for i in range (1, k+1):
		c = c*(n- ( k-i) )
		c = c // i
	return c 

print binomialcoefficient (20,10)


