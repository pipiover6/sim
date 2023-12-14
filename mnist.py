# Digit classifier [aka the 'Hello, World!' of Machine Learning]
# I'll follow the tutorial of "Mr. P Solver" at https://www.youtube.com/watch?v=gBw0u_5u0qU

# Download the mnist dataset at [www.di.ens.fr/~lelarge/MNIST.tar.gz] unzip it and save the the path
mnist_path = '/mnt/c/Users/pipiover6/Downloads/MNIST/MNIST/processed/'

import torch
import torch.nn as nn
from torch.optim import SGD
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
import torchvision
import numpy as np
import matplotlib.pyplot as plt

class DigitsDataset(Dataset):
    def __init__(self, filepath):
        self.x, self.y = torch.load(filepath)    # x is the vector of images of digits. Each image is 28x28 pixels. Each pixel is an integer in 0-255 [one byte]. 
                                                 # y is the vector of labels (actual digit) drawn the images.
        self.x = self.x / 255.                   # we normalize pixels to floats between 0 and 1.
        self.y = F.one_hot(self.y, num_classes=10).to(float)    # we normalize y so that instead of a digit (like 2) we keep an indicator vector of length 10 (like [0,0,1,0,0,...])
    
    # returns the number of (labeled) images in the dataset
    def __len__(self): 
        return self.x.shape[0]
      
    def __getitem__(self, ix): 
        return self.x[ix], self.y[ix]

# The mnist dataset contains 60k images for training and 10k for testing
train_ds = DigitsDataset(mnist_path+'training.pt')
assert(train_ds.x.shape == torch.Size([60000, 28, 28]))
test_ds = DigitsDataset(mnist_path+'test.pt')
assert(train_ds.y.shape == torch.Size([60000, 10]))

# we process images in batches (instead of one by one)
batch_size = 5
train_dl = DataLoader(train_ds, batch_size=batch_size)
assert(len(train_dl) * batch_size == len(train_ds))


class MyNeuralNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.Matrix1 = nn.Linear(28**2,100)  # first layer has 28*28 nodes (one for each pixel)
        self.Matrix2 = nn.Linear(100,50)     # second layer has 100
        self.Matrix3 = nn.Linear(50,10)      # third has 50 and last has 10 (one for each possible digit)
        self.R = nn.ReLU()                   # our model is so far linear - to give it tree like branching, we use a rectifier after each layer
    def forward(self,x):
        x = x.view(-1,28**2)  # this flattens a 2d image to a vector
        x = self.R(self.Matrix1(x))
        x = self.R(self.Matrix2(x))
        x = self.Matrix3(x)
        return x.squeeze()

f = MyNeuralNet()  # TA DA! Here we have our classifier. We only need to train it...

n_epochs = 10 # we relearn our model's weights going over the train data several times
def train_model(dl, f, n_epochs=n_epochs):
    # Optimization
    opt = SGD(f.parameters(), lr=0.01)  # we seek to minimize the loss using gradient descent
    # we measure our prediction [loss] using cross entropy (fit for classification)
    L = nn.CrossEntropyLoss()

    # Train model
    losses = []
    epochs = []
    for epoch in range(n_epochs):
        print(f'Epoch {epoch}')
        N = len(dl)
        for i, (x, y) in enumerate(dl):
            # Update the weights of the network
            opt.zero_grad() 
            loss_value = L(f(x), y) 
            loss_value.backward() 
            opt.step() 
            # Store training data
            epochs.append(epoch+i/N)
            losses.append(loss_value.item())
    print("Done training!")
    return np.array(epochs), np.array(losses)


epoch_data, loss_data = train_model(train_dl, f) # we trained our model! 

#we can graph the loss over the epochs
def plot_epoch_loss(epoch_data, loss_data):
  plt.plot(epoch_data, loss_data)
  plt.xlabel('Epoch Number')
  plt.ylabel('Cross Entropy')
  plt.title('Cross Entropy (per batch)')
  plt.show()

def model_correctness_percentage(on = "test"):
  ds = test_ds if (on == "test") else train_ds
  correct_count = 0
  n_samples = len(ds)
  incorrect_labeling_example = None
  for i in range(n_samples):
    prediction = f(ds.x[i]).argmax()
    target = ds.y[i].argmax()
    if target == prediction:
      correct_count += 1
    else:
        incorrect_labeling_example = (i, ds.x[i], target, prediction)
  return (correct_count / n_samples, incorrect_labeling_example)

train_prcn = model_correctness_percentage(on="train")[0]
print(f"Model correctness on training {train_prcn}")
prcn, bad = model_correctness_percentage()
print(f"Model correctness on test {prcn}")
if bad:
    i,img,target,prediction = bad
    print(f"Example of bad labeling: test data with index {i}, target = {target}, prediction = {prediction}")
    plt.imshow(img)
    plt.show()
