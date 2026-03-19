import pandas as pd
import torch
import torch.nn as nn
import pdb

# Define the autoencoder architecture
# Define the autoencoder architecture
# Define the autoencoder architecture



class Encoder(nn.Module):
    
    def __init__(self, encoded_space_dim):
        super().__init__()
        
        ### Convolutional section
        self.encoder_cnn = nn.Sequential(
            nn.Conv2d(3, 8, kernel_size = 3, stride=1, padding=1),
            nn.MaxPool2d(kernel_size=2),
            nn.ReLU(True),
            nn.Conv2d(8, 16, kernel_size = 3, stride=1, padding=1),
            nn.MaxPool2d(kernel_size=2),
            nn.ReLU(True),
            nn.Conv2d(16, 32, kernel_size = 3, stride=1, padding=1),
            nn.MaxPool2d(kernel_size=2),
            nn.ReLU(True),
            nn.Conv2d(32, 64, kernel_size = 3, stride=1, padding=1),
            nn.MaxPool2d(kernel_size=2),
            nn.ReLU(True),
            nn.Conv2d(64, 128, kernel_size = 3, stride=1, padding=1),
            nn.MaxPool2d(kernel_size=2),
            nn.ReLU(True)
        )
        
        ### Flatten layer
        self.flatten = nn.Flatten(start_dim=1)
### Linear section
        self.encoder_lin = nn.Sequential(
            nn.Linear(8 * 8 * 128, encoded_space_dim),
            # nn.ReLU(True),
            # nn.Linear(128, encoded_space_dim)
        )
        
    def forward(self, x):
        # pdb.set_trace()
        x = self.encoder_cnn(x)
        x = self.flatten(x)
        x = self.encoder_lin(x)
        return x


class Decoder(nn.Module):
    def __init__(self, encoded_space_dim):
        super().__init__()
        self.decoder_lin = nn.Sequential(
            nn.Linear(encoded_space_dim, 8 * 8 * 128),
            nn.ReLU(True),
            # nn.Linear(128, 8 * 8 * 128),
            # nn.ReLU(True)
        )

        self.unflatten = nn.Unflatten(dim=1, unflattened_size=(128, 8, 8))

        self.decoder_conv = nn.Sequential(
            nn.Upsample(scale_factor=2),
            nn.Conv2d(in_channels=128, out_channels=64, kernel_size=3, stride=1, padding=1),
            nn.ReLU(),
            nn.Upsample(scale_factor=2),
            nn.Conv2d(in_channels=64, out_channels=32, kernel_size=3, stride=1, padding=1),
            nn.ReLU(),
            nn.Upsample(scale_factor=2),
            nn.Conv2d(in_channels=32, out_channels=16, kernel_size=3, stride=1, padding=1),
            nn.ReLU(),
            nn.Upsample(scale_factor=2),
            nn.Conv2d(in_channels=16, out_channels=8, kernel_size=3, stride=1, padding=1),
            nn.ReLU(),
            nn.Upsample(scale_factor=2),
            nn.Conv2d(in_channels=8, out_channels=3, kernel_size=3, stride=1, padding=1),
            # nn.Sigmoid(),  # Use Sigmoid activation for pixel values between 0 and 1
            nn.Tanh(),     # I am using this activation because I am normalizing my image between -1 and 1
        )

    def forward(self, x):
        x = self.decoder_lin(x)
        x = self.unflatten(x)
        x = self.decoder_conv(x)
        return x
    
class  Autoencoder(nn.Module):
    def __init__(self, encoded_space_dim):
        super().__init__()
        self.encoder = Encoder(1024)
        self.decoder = Decoder(1024)
        
    def forward(self, x):
        bottelneck = self.encoder(x)
        decode = self.decoder(bottelneck)
        
        return decode, bottelneck



# Create an instance of the autoencoder

if __name__ == '__main__':
    pdb.set_trace()
    autoencoder = Autoencoder()
    input_rand = torch.rand(3, 256, 256)
    out, bottleneck = autoencoder(input_rand)


