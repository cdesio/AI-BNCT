import os
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

BNCT_PRIMARY_GEOMETRY_FEATURES = [
    "InitialEnergy_MeV",
    "SourceDistance_um",
    "Energy",
    "PosX_um",
    "PosY_um",
    "PosZ_um",
    "StopPosX_um",
    "StopPosY_um",
    "StopPosZ_um",
    "TraLen_cell_um",
    "TraLen_chro_um",
    "StartDistToOrigin_um",
    "StopDistToOrigin_um",
    "DeltaRadius_um",
    "TrackLength_um",
    "DirectionX",
    "DirectionY",
    "DirectionZ",
    "RadialDirectionCosine",
    "ClosestApproachT",
    "ClosestDistanceToOrigin_um",
    "HasDamage",
]


@dataclass
class Config:
    DATA_MODE: str = "csv_damage"
    DATA_GLOB: str = "/Users/yw18581/work/AI-BNCT/ForChiara-OldData/CSVfiles/*.csv"
    PRIMARY_SUMMARY_PATH: str = "data/bnct_primary_summary.csv"

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


@dataclass
class PrimaryGeometryConfig(Config):
    DATA_MODE: str = "primary_summary"
    PRIMARY_SUMMARY_PATH: str = "data/bnct_primary_summary.csv"
    TOTAL_EVENTS: int = 1000000
    FEATURE_NAMES: list[str] = field(default_factory=lambda: BNCT_PRIMARY_GEOMETRY_FEATURES.copy())
    BINARY_FEATURE_NAMES: list[str] = field(default_factory=lambda: ["HasDamage"])
    DAMAGE_FEATURE_NAME: str = "HasDamage"
    TRANSFORMER_LOC: str = "transformers/quantile_transformers_bnct_primary_geometry.pkl"
    RUN_TAG: str = "bnct_primary_geometry_steg"


def make_config():
    config_name = os.environ.get("BNCT_CONFIG", "csv_damage").strip().lower()
    if config_name in {"csv_damage", "csv", ""}:
        config = Config()
    elif config_name in {"primary_geometry", "primary", "geometry"}:
        config = PrimaryGeometryConfig()
    else:
        raise ValueError(
            "Unknown BNCT_CONFIG. Use 'csv_damage' or 'primary_geometry'. "
            f"Got {config_name!r}."
        )

    overrides = {
        "BNCT_DEVICE": ("DEVICE", str),
        "BNCT_EPOCHS": ("EPOCHS", int),
        "BNCT_TOTAL_EVENTS": ("TOTAL_EVENTS", int),
        "BNCT_CHECKPOINT_INTERVAL": ("CHECKPOINT_INTERVAL", int),
        "BNCT_NOISE_STEPS": ("NOISE_STEPS", int),
        "BNCT_LR": ("LR", float),
        "BNCT_PRIMARY_SUMMARY_PATH": ("PRIMARY_SUMMARY_PATH", str),
        "BNCT_RUN_TAG": ("RUN_TAG", str),
    }
    for env_name, (attribute, cast) in overrides.items():
        value = os.environ.get(env_name)
        if value is not None:
            setattr(config, attribute, cast(value))
    return config
