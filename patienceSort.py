#todo implement using heap/sorted tree instead of sorted list

def getMaxDeck(decks):
	assert(decks)
	curMaxInd = 0
	curMax = decks[0][-1]
	for i in range(len(decks)):
		assert(decks[i])
		if decks[i][-1] >= curMax:
			curMaxInd = i
			curMax = decks[i][-1]
	return (curMaxInd, curMax)
		

def patienceSort(arr):
	decks = []
	for a in arr:
		inserted = False
		for d in decks:
			assert(d)
			assert(not inserted)
			if a >= d[-1]:
				d.append(a)
				inserted = True
				break
		if (not inserted):
			decks.append([a])

	res = []
	for _ in range(len(arr)):
		print(decks)
		(i,a) = getMaxDeck(decks)
		res.append(a)
		decks[i].pop()
		if (not decks[i]):
			del decks[i]
	return res
