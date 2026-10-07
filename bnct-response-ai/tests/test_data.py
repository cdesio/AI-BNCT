import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

from bnct_response.data import _aggregate_steps, _assign_splits
from bnct_response.ps_two_step import _f1_threshold, _probability_logit
from bnct_response.steg import inverse_steg_output


def test_aggregate_steps_uses_last_post_energy_and_sums_deposition():
    steps = pd.DataFrame({
        "SeedID": [1, 1], "EventID": [2, 2], "TrackID": [1, 1],
        "VoxelID": [3, 3], "StepNumber": [2, 1],
        "KEPost_eV": [70.0, 90.0], "Edep_eV": [20.0, 10.0],
        "StepLength_nm": [2.0, 1.0],
    })
    result = _aggregate_steps(steps).iloc[0]
    assert result.ExitEnergy_eV == 70.0
    assert result.Edep_eV == 30.0
    assert result.TrackLength_nm == 3.0
    assert result.StepCount == 2
    assert np.isclose(result.MeanLET_keV_um, 10.0)


def test_splits_never_separate_one_primary():
    frame = pd.DataFrame([
        {"case_id": "a", "SeedID": 1, "EventID": event, "VoxelID": voxel}
        for event in range(20) for voxel in range(3)
    ])
    result = _assign_splits(frame)
    counts = result.groupby(["case_id", "SeedID", "EventID"])["split"].nunique()
    assert counts.eq(1).all()
    assert set(result.split) == {"train", "val", "test"}


def test_two_step_threshold_selects_best_f1_cutoff():
    truth = np.array([0, 0, 1, 1])
    probability = np.array([0.1, 0.4, 0.6, 0.9])
    threshold = _f1_threshold(truth, probability)
    assert 0.4 < threshold <= 0.6


def test_probability_logit_is_finite_at_boundaries():
    values = _probability_logit(np.array([0.0, 0.5, 1.0]))
    assert np.isfinite(values).all()
    assert values[0, 0] < 0 < values[2, 0]


def test_log_standard_steg_transform_round_trip():
    values = np.array([[0.0, 1.0], [3.0, 8.0]])
    transformer = StandardScaler().fit(np.log1p(values))
    scaled = transformer.transform(np.log1p(values))
    restored = inverse_steg_output(
        scaled, transformer, {"output_transform": "log_standard"}
    )
    assert np.allclose(restored, values)
