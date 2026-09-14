import torch
import numpy as np
import pickle


def generate_momenta(model, config, diffusion, quantile_transformers=None, batch_size=512):
    device = config.DEVICE

    # Generate scaled samples
    samples = diffusion.sample(model, batch_size)

    fake_scaled = samples.clone()

    # Convert to real-space if scaling exists
    samples_np = samples.cpu().numpy()

    if config.INPUT_SCALER_FLAG:
        for i, target in enumerate(quantile_transformers.keys()):
              samples_np[:,i] = quantile_transformers[target].inverse_transform(samples_np[:,i].reshape(-1,1)).flatten()

    binary_features = getattr(config, "BINARY_FEATURE_NAMES", [])
    for feature in binary_features:
        if feature in config.FEATURE_NAMES:
            feature_index = config.FEATURE_NAMES.index(feature)
            samples_np[:, feature_index] = (samples_np[:, feature_index] >= 0.5).astype(float)

    fake_unscaled = torch.tensor(samples_np, dtype=torch.float32, device=device)

    return fake_unscaled, fake_scaled
