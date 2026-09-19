"""圆形 DNA 双酶切核心算法。

纯 Python 实现，不依赖任何外部服务或第三方库：
- 按圆环相邻切点计算片段长度（含跨原点片段）；
- 按分辨容差把两种样本的片段聚为同一条带；
- 枚举酶对并按「可区分条带数 → 总条带数 → 酶名」排序。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from typing import Mapping, Sequence

MIN_CUTS_PER_SAMPLE = 2  # 环状 DNA 至少需要两个切点才能产生可比较的酶切片段


# ---------------------------------------------------------------------------
# 数据结构
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Fragment:
    """两个相邻切点之间的片段。跨原点片段的 end < start。"""

    start: int
    end: int
    length: int
    crosses_origin: bool
    visible: bool


@dataclass
class Band:
    """容差内聚成的一个条带，成员可来自参考型与变异型。"""

    min_length: float
    max_length: float
    reference_lengths: list[float] = field(default_factory=list)
    variant_lengths: list[float] = field(default_factory=list)
    center: float = 0.0

    @property
    def reference_count(self) -> int:
        return len(self.reference_lengths)

    @property
    def variant_count(self) -> int:
        return len(self.variant_lengths)

    @property
    def distinguishable(self) -> bool:
        """两种样本在该条带上的片段数不同 → 可用于鉴别。"""
        return self.reference_count != self.variant_count


@dataclass
class SampleDigest:
    """单一样本在某个酶对下的酶切结果。"""

    cut_positions: list[int]
    fragments: list[Fragment]
    band_count: int = 0


@dataclass
class PairResult:
    """一个酶对的完整评估结果。"""

    enzymes: tuple[str, str]
    reference: SampleDigest
    variant: SampleDigest
    bands: list[Band]
    distinguishable_bands: int
    total_bands: int

    @property
    def discriminates(self) -> bool:
        return self.distinguishable_bands > 0


@dataclass
class SkippedPair:
    enzymes: tuple[str, str]
    reason: str


@dataclass
class RecommendationResult:
    pairs: list[PairResult]
    skipped: list[SkippedPair]
    enzyme_pair_count: int


# ---------------------------------------------------------------------------
# 片段计算
# ---------------------------------------------------------------------------
def compute_fragments(
    plasmid_length: int,
    cut_positions: Sequence[int],
    min_visible_fragment: float = 0,
) -> list[Fragment]:
    """按圆环相邻切点计算片段长度。

    切点排序后，相邻切点之差即为片段长度；最后一个切点绕回第一个切点
    的片段为跨原点片段（crosses_origin=True）。
    """
    points = sorted(set(cut_positions))
    if len(points) < MIN_CUTS_PER_SAMPLE:
        raise ValueError(
            f"有效切点不足 {MIN_CUTS_PER_SAMPLE} 个，无法对环状 DNA 进行双酶切分析"
        )

    fragments: list[Fragment] = []
    for start, end in zip(points, points[1:]):
        length = end - start
        fragments.append(
            Fragment(
                start=start,
                end=end,
                length=length,
                crosses_origin=False,
                visible=length >= min_visible_fragment,
            )
        )
    wrap_length = plasmid_length - points[-1] + points[0]
    fragments.append(
        Fragment(
            start=points[-1],
            end=points[0],
            length=wrap_length,
            crosses_origin=True,
            visible=wrap_length >= min_visible_fragment,
        )
    )
    return fragments


# ---------------------------------------------------------------------------
# 条带聚类
# ---------------------------------------------------------------------------
def cluster_bands(
    reference_lengths: Sequence[float],
    variant_lengths: Sequence[float],
    tolerance: float,
) -> list[Band]:
    """把两种样本的可见片段按长度升序、容差内聚为同一条带。

    单链聚类：排序后相邻片段长度差 <= tolerance 即归入当前条带，
    否则开启新条带。条带中心取成员长度的算术平均。
    """
    tagged = [(float(x), "reference") for x in reference_lengths]
    tagged += [(float(x), "variant") for x in variant_lengths]
    tagged.sort(key=lambda item: item[0])

    bands: list[Band] = []
    for length, source in tagged:
        if bands and length - bands[-1].max_length <= tolerance:
            band = bands[-1]
            band.max_length = length
        else:
            band = Band(min_length=length, max_length=length)
            bands.append(band)
        if source == "reference":
            band.reference_lengths.append(length)
        else:
            band.variant_lengths.append(length)

    for band in bands:
        members = band.reference_lengths + band.variant_lengths
        band.center = sum(members) / len(members)
    return bands


# ---------------------------------------------------------------------------
# 酶对评估与推荐
# ---------------------------------------------------------------------------
def _merged_cuts(
    sites: Mapping[str, Sequence[int]], enzyme_a: str, enzyme_b: str
) -> list[int]:
    """合并两种酶在某样本中的切点（同一位置只切一次，取并集）。"""
    return sorted(set(sites.get(enzyme_a, ())) | set(sites.get(enzyme_b, ())))


def evaluate_pair(
    plasmid_length: int,
    enzyme_a: str,
    enzyme_b: str,
    reference_sites: Mapping[str, Sequence[int]],
    variant_sites: Mapping[str, Sequence[int]],
    min_visible_fragment: float,
    band_tolerance: float,
) -> PairResult | None:
    """评估一个酶对；任一样本有效切点不足两个时返回 None（拒绝该酶对）。"""
    ref_cuts = _merged_cuts(reference_sites, enzyme_a, enzyme_b)
    var_cuts = _merged_cuts(variant_sites, enzyme_a, enzyme_b)
    if len(ref_cuts) < MIN_CUTS_PER_SAMPLE or len(var_cuts) < MIN_CUTS_PER_SAMPLE:
        return None

    ref_fragments = compute_fragments(plasmid_length, ref_cuts, min_visible_fragment)
    var_fragments = compute_fragments(plasmid_length, var_cuts, min_visible_fragment)

    ref_visible = [f.length for f in ref_fragments if f.visible]
    var_visible = [f.length for f in var_fragments if f.visible]
    bands = cluster_bands(ref_visible, var_visible, band_tolerance)

    reference = SampleDigest(cut_positions=ref_cuts, fragments=ref_fragments)
    variant = SampleDigest(cut_positions=var_cuts, fragments=var_fragments)
    reference.band_count = sum(1 for b in bands if b.reference_count > 0)
    variant.band_count = sum(1 for b in bands if b.variant_count > 0)

    return PairResult(
        enzymes=(enzyme_a, enzyme_b),
        reference=reference,
        variant=variant,
        bands=bands,
        distinguishable_bands=sum(1 for b in bands if b.distinguishable),
        total_bands=len(bands),
    )


def recommend(
    plasmid_length: int,
    reference_sites: Mapping[str, Sequence[int]],
    variant_sites: Mapping[str, Sequence[int]],
    min_visible_fragment: float,
    band_tolerance: float,
) -> RecommendationResult:
    """枚举全部酶对，评估并按推荐优先级排序。

    排序规则：可区分条带数降序 → 总条带数升序（图谱越简单越好）→ 酶名升序。
    任一样本有效切点不足两个的酶对会被拒绝并记录在 skipped 中。
    """
    enzymes = sorted(set(reference_sites) | set(variant_sites))
    pairs: list[PairResult] = []
    skipped: list[SkippedPair] = []

    for enzyme_a, enzyme_b in combinations(enzymes, 2):
        result = evaluate_pair(
            plasmid_length,
            enzyme_a,
            enzyme_b,
            reference_sites,
            variant_sites,
            min_visible_fragment,
            band_tolerance,
        )
        if result is None:
            skipped.append(
                SkippedPair(
                    enzymes=(enzyme_a, enzyme_b),
                    reason=(
                        f"任一样本在该酶对下的有效切点不足 {MIN_CUTS_PER_SAMPLE} 个，已拒绝"
                    ),
                )
            )
        else:
            pairs.append(result)

    pairs.sort(key=lambda r: (-r.distinguishable_bands, r.total_bands, r.enzymes))
    return RecommendationResult(
        pairs=pairs, skipped=skipped, enzyme_pair_count=len(list(combinations(enzymes, 2)))
    )
