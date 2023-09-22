# We have a stick and break it into 3 parts by selecting two uniform division points on the stick. 
# What is the probability we can form a triangle from the 3 parts?

# The full stick is [0,1]
def validTriangle(x,y):
    return (min(x,y) < 0.5 < max(x,y)) and abs(x-y) < 0.5
  
def probValidTriMC(num_trials = 50000):
    s = 0
    for _ in range(num_trials):
        x = np.random.uniform()
        y = np.random.uniform()
        s += int(validTriangle(x,y))
    return s/num_trials
