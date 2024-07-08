# put probabilities in [0,0.5] and divide to two (approx) equal groups
# note that the 0.5 biased toss is independent noise
pr1 = [(25 + i)/100 for i in range(13)]
pr2 = [(49 - i)/100 for i in range(12)]


# return first index of arr with value >= x, None if all are < x
def bin_search(arr, x):
    lo = 0
    hi = len(arr) - 1
    while True:
        if lo > hi:
            assert lo == hi + 1
            return None
        if lo == hi:
            if arr[lo] < x:
                return None
            return lo        

        mid = (lo + hi) // 2
        assert lo <= mid <= hi
        if arr[mid] >= x:
            hi = mid
        else:
            lo = mid + 1


# returns an array of length 3 ^ the number of input probabilities, so in our example 3^12
# each entry has the form [p,q,p/q,k,sum p,sum q] to be explained in the code
def prepare_arr(probs):
    ret = []
    M = len(probs)
    for k in range(3**M):
        p = 1
        q = 1
        for i in range(M):
            x = (k // 3**i) % 3

            if x == 0:
                p *= (1-probs[i])**2
                q *= (probs[i])**2

            if x == 2:
                p *= (probs[i])**2
                q *= (1-probs[i])**2

            if x==1:
                p *= 2 * probs[i] * (1-probs[i])
                q *= 2 * probs[i] * (1-probs[i])

        ret.append([p, q, p/q, k])
    ret = sorted(ret, key = lambda z : z[2])
    # backwards running inclusive sum of p
    s = 0
    for i in range(len(ret)):
        elm = ret[len(ret) - 1 - i]
        s += elm[0]
        elm.append(s)

    # farwards running exclusive sum of q
    s = 0
    for i in range(len(ret)):
        elm = ret[i]
        elm.append(s)
        s += elm[1]

    return ret


def solve(probs1, probs2):
    arr1 = prepare_arr(probs1)
    arr2 = prepare_arr(probs2)
    quotients2 = [j[2] for j in arr2]

    s = 0
    for x in arr1:
        p = x[0]
        q = x[1]
        i = bin_search(quotients2, q/p)
        
        if i != None:
            s += p * arr2[i][4] + q * arr2[i][5]
        else:
            s += q
            
    return s/2

print(solve(pr1, pr2))
