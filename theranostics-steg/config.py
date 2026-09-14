from dataclasses import dataclass, field


THERANOSTICS_PHOTON_FEATURES = [
    "x_mm",
    "y_mm",
    "z_mm",
    "px",
    "py",
    "pz",
    "kineticEnergy_MeV",
    "time_s",
]


@dataclass
class Config:
    DATA_GLOB: str = (
        "/Users/yw18581/work/theranostics/"
        "rn222_repetition_runs/rn222_2mmfully_10reps_10k/outputs/*_photon_boundaries.csv"
    )

    BATCH_SIZE: int = None
    STEP_SIZE: int = 50000
    MAX_FILES: int = 50
    TOTAL_EVENTS: int = 100000

    EPOCHS: int = 30
    LR: float = 1e-5
    DEVICE: str = "cpu"

    @property
    def TRAIN_SPLIT(self):
        return int(0.7 * self.TOTAL_EVENTS)

    @property
    def VAL_SPLIT(self):
        return int(0.85 * self.TOTAL_EVENTS)

    CHECKPOINT_INTERVAL: int = 10

    FEATURE_NAMES: list[str] = field(default_factory=lambda: THERANOSTICS_PHOTON_FEATURES.copy())

    @property
    def FEATURES_DIM(self):
        return len(self.FEATURE_NAMES)

    TRANSFORMER_LOC: str = "transformers/quantile_transformers_theranostics_photons.pkl"
    INPUT_SCALER_FLAG: bool = True

    NOISE_STEPS: int = 100
    BETA_SCHEDULE: str = "linear"
    BETA_START: float = 1e-4
    BETA_END: float = 0.02

    MODEL_TYPE: str = "resnet"
    ATTENTION_GNN: bool = False
    HIDDEN_DIM: int = 256
    N_LAYERS: int = 5
    ACTIVATION: str = "relu"
    USE_LAYERNORM: bool = False
    TIME_CONDITIONING: str = "film"
    TIME_EMBED_DIM: int = 128

    RUNS_ROOT: str = "runs"
    RUN_TAG: str = "theranostics_test_runs"
    PLOT_2D_FEATURE_PAIRS: bool = False
