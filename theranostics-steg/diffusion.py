import torch
import math


class Diffusion:

    def __init__(self, config):
        self.device = config.DEVICE
        self.n_steps = config.NOISE_STEPS
        self.features_dim = config.FEATURES_DIM

        if config.BETA_SCHEDULE == "linear":
            self.beta = torch.linspace(
                config.BETA_START,
                config.BETA_END,
                self.n_steps,
                device=self.device
            )

        elif config.BETA_SCHEDULE == "cosine":
            self.beta = self.cosine_beta_schedule(self.n_steps).to(self.device)

        else:
            raise ValueError("Unknown beta schedule")

        self.alpha = 1.0 - self.beta
        self.alpha_hat = torch.cumprod(self.alpha, dim=0)

    def cosine_beta_schedule(self, timesteps, s=0.008):
        steps = timesteps + 1
        x = torch.linspace(0, timesteps, steps)
        alphas_cumprod = torch.cos(
            ((x / timesteps) + s) / (1 + s) * math.pi * 0.5
        ) ** 2

        alphas_cumprod = alphas_cumprod / alphas_cumprod[0]

        betas = 1 - (alphas_cumprod[1:] / alphas_cumprod[:-1])
        return torch.clamp(betas, 1e-5, 0.999)

    def sample_timesteps(self, batch_size):
        return torch.randint(0, self.n_steps, (batch_size,), device=self.device)

    def noise_data(self, x0, t):
        sqrt_alpha_hat = torch.sqrt(self.alpha_hat[t])[:, None]
        sqrt_one_minus_alpha_hat = torch.sqrt(1 - self.alpha_hat[t])[:, None]

        noise = torch.randn_like(x0)
        xt = sqrt_alpha_hat * x0 + sqrt_one_minus_alpha_hat * noise

        return xt, noise

    def training_step(self, model, x0):

        t = self.sample_timesteps(x0.shape[0])
        xt, noise = self.noise_data(x0, t)

        t_norm = t.float() / self.n_steps

        pred_noise = model(xt, t_norm)

        loss = ((noise - pred_noise) ** 2).mean()

        return loss
    
    def sample(self, model, n_samples):
        model.eval()
        x = torch.randn((n_samples, self.features_dim), device=self.device)

        for i in reversed(range(self.n_steps)):
            t = torch.ones(n_samples, device=self.device) * (i / self.n_steps)

            beta = self.beta[i]
            alpha = self.alpha[i]
            alpha_hat = self.alpha_hat[i]

            pred_noise = model(x, t)

            # DDPM reverse step
            x = (1 / torch.sqrt(alpha)) * (
                x - (beta / torch.sqrt(1 - alpha_hat)) * pred_noise
            )

            # Add noise except at final step
            if i > 0:
                noise = torch.randn_like(x)
                x += torch.sqrt(beta) * noise

        return x