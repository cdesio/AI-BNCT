from dataclasses import dataclass, field


BNCT_TRACK_FEATURES = [
    "PosX_um",
    "PosY_um",
    "PosZ_um",
    "StopPosX_um",
    "StopPosY_um",
    "StopPosZ_um",
    "Dist_Start_to_Origin",
    "Dist_Stop_to_Start",
]

BNCT_RAW_DAMAGE_COLUMNS = [
    "SSB",
    "DSB",
    "DSBp",
    "DSBpp",
]

BNCT_DAMAGE_INDICATOR_FEATURES = [
    "Damage",
]

BNCT_DAMAGE_FEATURES = [
    *BNCT_TRACK_FEATURES,
    *BNCT_DAMAGE_INDICATOR_FEATURES,
]


@dataclass
class Config:
    DATA_GLOB: str = "/Users/yw18581/work/AI-BNCT/ForChiara-OldData/CSVfiles/*.csv"

    BATCH_SIZE: int = None
    STEP_SIZE: int = 50000
    MAX_FILES: int = None

    EPOCHS: int = 100
    LR: float = 1e-5

    DEVICE: str = "cpu"
    TOTAL_EVENTS: int = 100000

    @property
    def TRAIN_SPLIT(self):
        return int(0.7 * self.TOTAL_EVENTS)

    @property
    def VAL_SPLIT(self):
        return int(0.85 * self.TOTAL_EVENTS)

    CHECKPOINT_INTERVAL: int = 10

    FEATURE_NAMES: list[str] = field(default_factory=lambda: BNCT_DAMAGE_FEATURES.copy())
    BINARY_FEATURE_NAMES: list[str] = field(
        default_factory=lambda: BNCT_DAMAGE_INDICATOR_FEATURES.copy()
    )
    DAMAGE_FEATURE_NAME: str = "Damage"

    @property
    def FEATURES_DIM(self):
        return len(self.FEATURE_NAMES)

    TRANSFORMER_LOC: str = "transformers/quantile_transformers_bnct_csv_damage_indicators.pkl"
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
    RUN_TAG: str = "bnct_csv_damage_smoke"
    PLOT_2D_FEATURE_PAIRS: bool = False
