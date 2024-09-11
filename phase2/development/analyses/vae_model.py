# Gemini generated code
import torch
from torch import nn
from torch.nn import functional as F

class VAE(nn.Module):
  def __init__(self, input_dim, latent_dim):
    super(VAE, self).__init__()
    self.encoder = nn.Sequential(
        nn.Linear(input_dim, 256),
        nn.ReLU(),
        nn.Linear(256, 128),
        nn.ReLU()
    )

    self.fc1 = nn.Linear(128, latent_dim)
    self.fc2 = nn.Linear(128, latent_dim)

    self.decoder = nn.Sequential(
        nn.Linear(latent_dim, 128),
        nn.ReLU(),
        nn.Linear(128, 256),
        nn.ReLU(),
        nn.Linear(256, input_dim)
    )

  def reparameterize(self, mu, logvar):
    std = torch.exp(0.5*logvar)
    eps = torch.randn_like(std)
    return mu + eps * std

  def forward(self, x):
    h = self.encoder(x)
    mu = self.fc1(h)
    logvar = self.fc2(h)
    z = self.reparameterize(mu, logvar)
    recon_x = self.decoder(z)
    return recon_x, mu, logvar

# Example Usage
vae = VAE(10, 4) # 10 input features, 4 latent dimensions
optimizer = torch.optim.Adam(vae.parameters(), lr=0.001)

# Train loop (replace with your data loader)
for epoch in range(10):
  for data in data_loader:
    # reconstruction loss
    recon_batch, mu, logvar = vae(data)
    bce_loss = F.binary_cross_entropy(recon_batch, data, reduction='sum')

    # KL divergence loss
    kld_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())

    # total loss
    loss = bce_loss + kld_loss
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

# Generate new data from latent space (optional)
with torch.no_grad():
  z = torch.randn(1, 4) # Sample random noise
  new_data = vae.decoder(z)
  print(new_data)


# also take a look at https://github.com/adrianjav/heterogeneous_vaes

# Perplexity generated VAE
import torch
import torch.nn as nn
from torch.distributions import Normal, Independent

class TabularVAE(nn.Module):
    def __init__(self, input_dim, latent_dim):
        super(TabularVAE, self).__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 512),
            nn.ReLU(),
            nn.Linear(512, 256),
            nn.ReLU(),
        )
        self.mu = nn.Linear(256, latent_dim)
        self.log_var = nn.Linear(256, latent_dim)
        
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 256),
            nn.ReLU(),
            nn.Linear(256, 512),
            nn.ReLU(),
            nn.Linear(512, input_dim),
        )

    def encode(self, x):
        h = self.encoder(x)
        return self.mu(h), self.log_var(h)

    def reparameterize(self, mu, log_var):
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z):
        return self.decoder(z)

    def forward(self, x):
        mu, log_var = self.encode(x)
        z = self.reparameterize(mu, log_var)
        x_recon = self.decode(z)
        return x_recon, mu, log_var