"""Shared feature and target definitions."""

CONDITION_COLUMNS = [
    "Particle",
    "InitialEnergy_MeV",
    "Distance_um",
    "InitialDirX",
    "InitialDirY",
    "InitialDirZ",
]

TRANSPORT_TARGETS = [
    "EntryEnergy_MeV",
    "ExitEnergy_MeV",
    "Edep_keV",
    "TrackLength_nm",
    "StepCount",
    "MeanLET_keV_um",
    "EntryDirX",
    "EntryDirY",
    "EntryDirZ",
    "LocalEntryX_nm",
    "LocalEntryY_nm",
    "LocalEntryZ_nm",
    "TransverseDisplacement_um",
]

DAMAGE_TRANSPORT_COLUMNS = [
    "EntryEnergy_MeV",
    "ExitEnergy_MeV",
    "Edep_keV",
    "TrackLength_nm",
    "StepCount",
    "MeanLET_keV_um",
]

# Damage is marginalized over unresolved microscopic entry position. Feeding a
# point prediction of the multimodal voxel-face coordinates creates an
# unphysical average entry state in the chained baseline.
DAMAGE_INPUT_COLUMNS = CONDITION_COLUMNS + DAMAGE_TRANSPORT_COLUMNS

DAMAGE_COUNT_TARGETS = [
    "TotalSB",
    "TotalSSB",
    "TotalCSSB",
    "TotalDSB",
    "TotalCDSB",
    "DirectSB",
    "IndirectSB",
]

DAMAGE_TARGETS = ["AnyDamage", "AnyDSB"] + DAMAGE_COUNT_TARGETS

PRIMARY_KEY = ["case_id", "SeedID", "EventID"]
ENCOUNTER_KEY = PRIMARY_KEY + ["TrackID", "VoxelID"]

# Local phase-space state available before Geant4-DNA replay. World position,
# voxel ID and primary distance are provenance, not response-model inputs.
PS_INPUT_COLUMNS = [
    "Particle",
    "EntryEnergy_MeV",
    "EntryDirX",
    "EntryDirY",
    "EntryDirZ",
    "LocalEntryX_nm",
    "LocalEntryY_nm",
    "LocalEntryZ_nm",
]

PS_HISTORY_INPUT_COLUMNS = PS_INPUT_COLUMNS + [
    "WorldEntryX_nm",
    "WorldEntryY_nm",
    "WorldEntryZ_nm",
    "Distance_um",
    "InitialEnergy_MeV",
    "InitialDirX",
    "InitialDirY",
    "InitialDirZ",
]

PS_RESPONSE_CONTINUOUS_TARGETS = ["DNAEdep_keV"]
PS_RESPONSE_COUNT_TARGETS = DAMAGE_COUNT_TARGETS
PS_GENERATIVE_COUNT_TARGETS = [
    "TotalSB", "TotalSSB", "TotalCSSB", "TotalDSB", "TotalCDSB"
]
PS_RESPONSE_TARGETS = PS_RESPONSE_CONTINUOUS_TARGETS + PS_GENERATIVE_COUNT_TARGETS
