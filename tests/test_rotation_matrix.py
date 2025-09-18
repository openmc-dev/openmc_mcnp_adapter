import numpy as np
import pytest

from openmc_mcnp_adapter import rotation_matrix


def is_orthogonal(R: np.ndarray, atol: float = 1e-12) -> bool:
    """Check if the matrix R is orthogonal"""
    I = np.identity(3)
    return np.allclose(R.T @ R, I, atol=atol) and np.allclose(R @ R.T, I, atol=atol)


@pytest.mark.parametrize(
    "v1, v2",
    [
        # Same direction, different magnitudes
        (np.array([1.0, 2.0, 3.0]), np.array([2.0, 4.0, 6.0])),
        (np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, 10.0])),
    ],
)
def test_rotation_parallel(v1, v2):
    """Test rotation_matrix for parallel vectors"""
    R = rotation_matrix(v1, v2)
    # Should be exactly identity from the parallel branch
    assert np.allclose(R, np.identity(3))
    assert is_orthogonal(R)
    assert np.isclose(np.linalg.det(R), 1.0)


@pytest.mark.parametrize(
    "v1, v2",
    [
        (np.array([0.0, 0.0, 1.0]), np.array([0.0, 0.0, -5.0])),  # z -> -z
        (np.array([1.0, 0.0, 0.0]), np.array([-2.0, 0.0, 0.0])),  # x -> -x
    ],
)
def test_rotation_antiparallel(v1, v2):
    """Test rotation_matrix for anti-parallel vectors"""
    R = rotation_matrix(v1, v2)

    # Maps v1 direction to v2 direction
    v2_hat = v2 / np.linalg.norm(v2)
    mapped = R @ (v1 / np.linalg.norm(v1))
    assert np.allclose(mapped, v2_hat, atol=1e-12)

    # No NaNs, orthogonal and det +1
    assert not np.isnan(R).any()
    assert is_orthogonal(R)
    assert np.isclose(np.linalg.det(R), 1.0, atol=1e-12)
    # 180-degree rotation has trace -1
    assert np.isclose(np.trace(R), -1.0, atol=1e-12)


@pytest.mark.parametrize(
    "v1, v2",
    [
        (np.array([0.0, 0.0, 1.0]), np.array([1.0, 2.0, 3.0])),
        (np.array([0.0, 1.0, 0.0]), np.array([3.0, -1.0, 2.0])),
    ],
)
def test_rotation_general(v1, v2):
    """Test rotation_matrix for general vectors"""
    R = rotation_matrix(v1, v2)

    # Maps v1 to v2 direction
    v2_hat = v2 / np.linalg.norm(v2)
    mapped = R @ (v1 / np.linalg.norm(v1))
    assert np.allclose(mapped, v2_hat, atol=1e-12)

    # Orthogonal and proper rotation
    assert is_orthogonal(R)
    assert np.isclose(np.linalg.det(R), 1.0, atol=1e-12)

    # A vector perpendicular to v1 remains perpendicular to R v1 (i.e., to v2_hat)
    # Use a simple perpendicular vector: pick the least-aligned basis axis and project out
    basis = [
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
        np.array([0.0, 0.0, 1.0])
    ]
    u1_hat = v1 / np.linalg.norm(v1)
    e = min(basis, key=lambda b: abs(np.dot(b, u1_hat)))
    x = e - np.dot(e, u1_hat) * u1_hat
    x /= np.linalg.norm(x)

    # Should be perpendicular to v2_hat
    assert np.isclose(np.dot(R @ x, v2_hat), 0.0, atol=1e-12)
