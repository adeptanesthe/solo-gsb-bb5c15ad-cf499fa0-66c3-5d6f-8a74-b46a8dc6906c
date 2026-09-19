"""核心算法单元测试：片段计算与条带聚类。"""
import pytest

from app.digest import cluster_bands, compute_fragments


def test_fragments_circular_adjacent_cuts():
    fragments = compute_fragments(6000, [3000, 500, 1200])
    lengths = [f.length for f in fragments]
    assert lengths == [700, 1800, 3500]
    assert sum(lengths) == 6000


def test_cross_origin_fragment_flagged():
    fragments = compute_fragments(6000, [500, 3000])
    cross = [f for f in fragments if f.crosses_origin]
    normal = [f for f in fragments if not f.crosses_origin]
    assert len(cross) == 1
    assert cross[0].start == 3000 and cross[0].end == 500
    assert cross[0].length == 6000 - 3000 + 500
    assert normal[0].length == 2500


def test_duplicate_cut_positions_merged():
    fragments = compute_fragments(6000, [500, 500, 3000])
    assert [f.length for f in fragments] == [2500, 3500]


def test_too_few_cuts_rejected():
    with pytest.raises(ValueError):
        compute_fragments(6000, [500])
    with pytest.raises(ValueError):
        compute_fragments(6000, [])


def test_min_visible_fragment_marks_invisible():
    fragments = compute_fragments(6000, [500, 900, 3000], min_visible_fragment=500)
    visibility = {f.length: f.visible for f in fragments}
    assert visibility == {400: False, 2100: True, 3500: True}


def test_cluster_bands_within_tolerance():
    # 1000/1005 差值在容差内 → 同一条带；3000 与 3100 超出容差 → 各自成带
    bands = cluster_bands([1000, 3000], [1005, 3100], tolerance=10)
    assert len(bands) == 3
    merged, ref_only, var_only = bands
    assert merged.reference_count == 1 and merged.variant_count == 1
    assert not merged.distinguishable
    assert ref_only.reference_count == 1 and ref_only.variant_count == 0
    assert ref_only.distinguishable
    assert var_only.reference_count == 0 and var_only.variant_count == 1
    assert var_only.distinguishable


def test_cluster_bands_single_linkage_chain():
    # 相邻差都在容差内 → 链式归入同一条带
    bands = cluster_bands([1000], [1008, 1015], tolerance=10)
    assert len(bands) == 1
    assert bands[0].reference_count == 1
    assert bands[0].variant_count == 2


def test_cluster_bands_gap_opens_new_band():
    bands = cluster_bands([1000], [1030], tolerance=10)
    assert len(bands) == 2
