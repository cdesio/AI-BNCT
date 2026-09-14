import torch
import torch.nn as nn
import torch.nn.functional as F


# -------------------------
# Utilities
# -------------------------

def get_activation(name):
    if name == "relu":
        return nn.ReLU()
    elif name == "silu":
        return nn.SiLU()
    elif name == "gelu":
        return nn.GELU()
    else:
        raise ValueError(f"Unknown activation: {name}")


# -------------------------
# Time Embedding
# -------------------------

class TimeEmbedding(nn.Module):
    def __init__(self, hidden_dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(1, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim)
        )

    def forward(self, t):
        t = t.unsqueeze(-1)
        return self.net(t)


# -------------------------
# MLP Block
# -------------------------

class MLPBlock(nn.Module):
    def __init__(self, dim, activation, use_layernorm):
        super().__init__()
        self.use_ln = use_layernorm
        if use_layernorm:
            self.ln = nn.LayerNorm(dim)

        self.fc1 = nn.Linear(dim, dim)
        self.fc2 = nn.Linear(dim, dim)
        self.act = get_activation(activation)

    def forward(self, x):
        h = x
        if self.use_ln:
            h = self.ln(h)

        h = self.fc1(h)
        h = self.act(h)
        h = self.fc2(h)
        return h


# -------------------------
# ResNet Block (supports time injection + FiLM)
# -------------------------

class ResBlock(nn.Module):
    def __init__(self, dim, config):
        super().__init__()

        self.use_ln = config.USE_LAYERNORM
        self.time_cond = config.TIME_CONDITIONING
        self.activation = get_activation(config.ACTIVATION)

        if self.use_ln:
            self.ln = nn.LayerNorm(dim)

        self.fc1 = nn.Linear(dim, dim)
        self.fc2 = nn.Linear(dim, dim)

        # FiLM conditioning
        if self.time_cond == "film":
            self.gamma = nn.Linear(config.TIME_EMBED_DIM, dim)
            self.beta = nn.Linear(config.TIME_EMBED_DIM, dim)

        if self.time_cond in ["per_block", "film"]:
            self.time_proj = nn.Linear(config.TIME_EMBED_DIM, dim)

    def forward(self, x, t_embed=None):
        h = x

        if self.use_ln:
            h = self.ln(h)

        # Time conditioning
        if self.time_cond == "per_block" and t_embed is not None:
            h = h + self.time_proj(t_embed)

        elif self.time_cond == "film" and t_embed is not None:
            gamma = self.gamma(t_embed)
            beta = self.beta(t_embed)
            h = gamma * h + beta

        h = self.fc1(h)
        h = self.activation(h)
        h = self.fc2(h)

        return x + h  # residual connection


# -------------------------
# GNN Layer
# -------------------------

class GNNLayer(nn.Module):
    def __init__(self, dim, config):
        super().__init__()

        self.use_ln = config.USE_LAYERNORM
        self.activation = get_activation(config.ACTIVATION)
        self.time_cond = config.TIME_CONDITIONING
        self.attention = config.ATTENTION_GNN

        if self.use_ln:
            self.ln = nn.LayerNorm(dim)

        self.msg = nn.Linear(dim, dim)
        self.update = nn.Linear(dim, dim)

        if self.time_cond == "film":
            self.gamma = nn.Linear(dim, dim)
            self.beta = nn.Linear(dim, dim)

        self.time_proj = nn.Linear(config.TIME_EMBED_DIM, dim)

    def forward(self, h, g=None, t_embed=None):
        if self.use_ln:
            h = self.ln(h)

        if self.time_cond == "per_block" and t_embed is not None:
            h = h + self.time_proj(t_embed).unsqueeze(1)

        elif self.time_cond == "film" and t_embed is not None:
            t_proj = self.time_proj(t_embed)
            gamma = self.gamma(t_proj).unsqueeze(1)
            beta = self.beta(t_proj).unsqueeze(1)
            h = gamma * h + beta

        # -------------------------
        # Global → node
        # -------------------------
        if g is not None:
            h = h + g.unsqueeze(1)

        if self.attention == False:
            # Treat all particles equally
            messages = self.msg(h)
            agg = messages.mean(dim=1, keepdim=True)

        else:
            # Attention GNN, learns weights for each particle
            messages = self.msg(h)  # (B, N, dim)

            # compute attention scores
            scores = torch.matmul(messages, messages.transpose(1,2)) / (messages.shape[-1] ** 0.5)
            weights = torch.softmax(scores, dim=-1)

            agg = torch.matmul(weights, messages).mean(dim=1, keepdim=True)

        h = h + self.update(agg)
        h = self.activation(h)

        if g is not None:
            g = g + agg.squeeze(1)

        return h, g


# -------------------------
# Main Model
# -------------------------

class MomentaDiffusionModel(nn.Module):
    def __init__(self, config):
        super().__init__()

        self.config = config
        self.model_type = config.MODEL_TYPE
        self.input_dim = config.FEATURES_DIM
        self.hidden = config.HIDDEN_DIM

        # Time embedding
        if config.TIME_CONDITIONING != "none":
            self.time_embed = TimeEmbedding(config.TIME_EMBED_DIM)
        else:
            self.time_embed = None

        # -------------------------
        # MLP / ResNet
        # -------------------------
        if self.model_type in ["mlp", "resnet"]:

            in_dim = self.input_dim
            if config.TIME_CONDITIONING == "input":
                in_dim += config.TIME_EMBED_DIM

            self.input_layer = nn.Linear(in_dim, self.hidden)

            blocks = []
            for _ in range(config.N_LAYERS):
                if self.model_type == "mlp":
                    blocks.append(
                        MLPBlock(self.hidden, config.ACTIVATION, config.USE_LAYERNORM)
                    )
                else:
                    blocks.append(ResBlock(self.hidden, config))

            self.blocks = nn.ModuleList(blocks)
            self.output_layer = nn.Linear(self.hidden, self.input_dim)

        # -------------------------
        # GNN
        # -------------------------
        elif self.model_type == "gnn":

            self.n_nodes = 4
            self.node_dim = 5   # px, py, pz, E, BPVIP
            self.global_dim = 3 # B vertex

            self.node_input = nn.Linear(self.node_dim, self.hidden)
            self.global_input = nn.Linear(self.global_dim, self.hidden)

            self.layers = nn.ModuleList([
                GNNLayer(self.hidden, config)
                for _ in range(config.N_LAYERS)
            ])

            self.node_output = nn.Linear(self.hidden, self.node_dim)
            self.global_output = nn.Linear(self.hidden, self.global_dim)

        else:
            raise ValueError(f"Unknown model type: {self.model_type}")

    # -------------------------
    # Forward
    # -------------------------
    def forward(self, x, t):

        # Time embedding
        if self.time_embed is not None:
            t_embed = self.time_embed(t)
        else:
            t_embed = None

        # -------------------------
        # MLP / ResNet
        # -------------------------
        if self.model_type in ["mlp", "resnet"]:

            if self.config.TIME_CONDITIONING == "input" and t_embed is not None:
                x = torch.cat([x, t_embed], dim=1)

            h = self.input_layer(x)

            for block in self.blocks:
                if self.model_type == "resnet":
                    h = block(h, t_embed)
                else:
                    h = h + block(h)  # optional residual-style MLP

            return self.output_layer(h)

        # -------------------------
        # GNN
        # -------------------------
        elif self.model_type == "gnn":

            B = x.shape[0]

            # -------------------------
            # Split input
            # -------------------------
            momenta = x[:, :16].view(B, 4, 4)
            bpvip   = x[:, 16:20].view(B, 4, 1)

            h = torch.cat([momenta, bpvip], dim=-1)  # (B, 4, 5)
            g = x[:, 20:23]                          # (B, 3)

            # -------------------------
            # Input projection
            # -------------------------
            h = self.node_input(h)
            g = self.global_input(g)

            # -------------------------
            # GNN layers
            # -------------------------
            for layer in self.layers:
                h, g = layer(h, g, t_embed)

            # -------------------------
            # Output
            # -------------------------
            h = self.node_output(h)
            g = self.global_output(g)

            return torch.cat([h.view(B, -1), g], dim=1)
