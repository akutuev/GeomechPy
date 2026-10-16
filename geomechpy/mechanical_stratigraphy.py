"""Mechanical stratigraphy: lithology classification from the Gamma Ray log.

Each depth sample is classified into a discrete lithology code from the Gamma
Ray (GR) reading using a single, user-defined gAPI cutoff (baseline): readings
at or above the cutoff are shale, otherwise the clean rock is treated as
sandstone. Coal and limestone cannot be resolved from GR alone, so they are
supplied as external flags that override the GR-based classification.

Lithology codes (see :class:`Lithology`): Sand = 0, Shale = 1, Limestone = 2,
Coal = 6.
"""

from enum import IntEnum


class Lithology(IntEnum):
    """Lithology codes used for mechanical stratigraphy classification.

    The integer values are the codes emitted by the lithology-determination
    functions and are intended to be written directly to a discrete log curve.
    """

    SAND = 0
    SHALE = 1
    LIMESTONE = 2
    COAL = 6


def determine_lithology(
    gamma_ray: float,
    gr_threshold: float = 75.0,
    coal_flag: bool = False,
    limestone_flag: bool = False,
) -> int:
    """Determine the mechanical stratigraphy lithology from a Gamma Ray reading.

    The Gamma Ray log is used to discriminate shale from non-shale rock: a reading
    at or above the threshold is classified as shale, otherwise the interval is
    treated as clean (sandstone). Coal and limestone are not resolvable from Gamma
    Ray alone, so they are supplied as external flags that override the Gamma Ray
    classification. The coal flag takes precedence over the limestone flag.

    The returned integer follows the mechanical stratigraphy convention encoded in
    :class:`Lithology`: Sand = 0, Shale = 1, Limestone = 2, Coal = 6.

    Args:
        gamma_ray (float): Gamma Ray reading. Unit: [gAPI]
        gr_threshold (float): Gamma Ray cutoff (baseline) separating clean (sand)
            from shale. Readings at or above this value are shale. Unit: [gAPI].
            Default: 75.0
        coal_flag (bool): When True the interval is classified as coal, overriding
            the Gamma Ray classification. Default: False
        limestone_flag (bool): When True the interval is classified as limestone,
            overriding the Gamma Ray classification. Ignored if ``coal_flag`` is
            True. Default: False

    Returns:
        int: Lithology code. Sand = 0, Shale = 1, Limestone = 2, Coal = 6
    """
    if coal_flag:
        return int(Lithology.COAL)
    if limestone_flag:
        return int(Lithology.LIMESTONE)
    if gamma_ray >= gr_threshold:
        return int(Lithology.SHALE)
    return int(Lithology.SAND)


def determine_lithology_array(
    gamma_ray: list[float],
    gr_threshold: float = 75.0,
    coal_flag: list[bool] | None = None,
    limestone_flag: list[bool] | None = None,
) -> list[int]:
    """Determine the mechanical stratigraphy lithology for a Gamma Ray log.

    Applies :func:`determine_lithology` to every sample of a Gamma Ray curve. The
    ``coal_flag`` and ``limestone_flag`` sequences, when provided, must be the same
    length as ``gamma_ray``; when omitted the corresponding flag is treated as
    False for every sample.

    Args:
        gamma_ray (list[float]): Gamma Ray readings. Unit: [gAPI]
        gr_threshold (float): Gamma Ray cutoff (baseline) separating clean (sand)
            from shale. Unit: [gAPI]. Default: 75.0
        coal_flag (list[bool] | None): Per-sample coal flags. Default: None (all False)
        limestone_flag (list[bool] | None): Per-sample limestone flags. Default: None (all False)

    Returns:
        list[int]: Lithology codes. Sand = 0, Shale = 1, Limestone = 2, Coal = 6

    Raises:
        ValueError: If ``coal_flag`` or ``limestone_flag`` is provided with a length
            that does not match ``gamma_ray``.
    """
    coal = coal_flag if coal_flag is not None else [False] * len(gamma_ray)
    limestone = limestone_flag if limestone_flag is not None else [False] * len(gamma_ray)

    if len(coal) != len(gamma_ray):
        raise ValueError("coal_flag must have the same length as gamma_ray")
    if len(limestone) != len(gamma_ray):
        raise ValueError("limestone_flag must have the same length as gamma_ray")

    return [
        determine_lithology(
            gamma_ray=gr,
            gr_threshold=gr_threshold,
            coal_flag=coal_sample,
            limestone_flag=limestone_sample,
        )
        for gr, coal_sample, limestone_sample in zip(gamma_ray, coal, limestone)
    ]
