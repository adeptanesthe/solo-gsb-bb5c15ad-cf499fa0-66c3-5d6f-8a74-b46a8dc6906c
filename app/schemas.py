"""请求 / 响应模型与领域校验。

位点越界、重复位点、酶数量不足等错误在此抛出，
由 FastAPI 统一转换为 422 响应。
"""
from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


# ---------------------------------------------------------------------------
# 请求
# ---------------------------------------------------------------------------
class DigestRequest(BaseModel):
    """双酶切鉴别请求。"""

    model_config = ConfigDict(str_strip_whitespace=True)

    plasmid_length: int = Field(..., gt=0, description="环状质粒总长度 (bp)")
    reference_sites: dict[str, list[int]] = Field(
        ..., description="参考型候选酶切位点：{酶名: [切点位置...]}"
    )
    variant_sites: dict[str, list[int]] = Field(
        ..., description="变异型候选酶切位点：{酶名: [切点位置...]}"
    )
    min_visible_fragment: float = Field(
        0, ge=0, description="最小可见片段长度 (bp)，短于该值的片段不计入条带"
    )
    band_tolerance: float = Field(
        0, ge=0, description="条带分辨容差 (bp)，长度差不超过该值的片段聚为同一条带"
    )
    top_n: int | None = Field(
        None, ge=1, le=100, description="最多返回的推荐方案数；缺省返回全部有效酶对"
    )

    @field_validator("reference_sites", "variant_sites", mode="before")
    @classmethod
    def _check_sites(cls, value):
        if not isinstance(value, dict):
            raise ValueError("酶切位点必须是 {酶名: [切点...]} 的字典")
        for enzyme, positions in value.items():
            name = str(enzyme).strip()
            if not name:
                raise ValueError("酶名不能为空")
            if not isinstance(positions, list):
                raise ValueError(f"酶 {name} 的切点必须是列表")
            seen: set[int] = set()
            for pos in positions:
                if isinstance(pos, bool) or not isinstance(pos, int):
                    raise ValueError(f"酶 {name} 的切点 {pos!r} 不是整数")
                if pos in seen:
                    raise ValueError(f"酶 {name} 存在重复位点 {pos}")
                seen.add(pos)
        return value

    @model_validator(mode="after")
    def _check_bounds_and_enzymes(self) -> "DigestRequest":
        enzymes = set(self.reference_sites) | set(self.variant_sites)
        if len(enzymes) < 2:
            raise ValueError("至少需要两种候选限制酶才能枚举酶对")
        for label, sites in (("参考型", self.reference_sites), ("变异型", self.variant_sites)):
            for enzyme, positions in sites.items():
                for pos in positions:
                    if not 0 <= pos < self.plasmid_length:
                        raise ValueError(
                            f"{label}样本酶 {enzyme} 的位点 {pos} 越界："
                            f"有效范围为 [0, {self.plasmid_length})"
                        )
        return self


# ---------------------------------------------------------------------------
# 响应
# ---------------------------------------------------------------------------
class FragmentOut(BaseModel):
    start: int = Field(..., description="片段起始切点")
    end: int = Field(..., description="片段终止切点；跨原点片段 end < start")
    length: int = Field(..., description="片段长度 (bp)")
    crosses_origin: bool = Field(..., description="是否为跨原点片段")
    visible: bool = Field(..., description="是否不小于最小可见片段")


class SampleDigestOut(BaseModel):
    cut_positions: list[int] = Field(..., description="合并两种酶后的有效切点（升序）")
    fragment_count: int
    fragments: list[FragmentOut]
    cross_origin_fragment: FragmentOut = Field(..., description="跨原点片段")
    visible_fragment_lengths: list[int] = Field(..., description="参与聚类的可见片段长度")
    band_count: int = Field(..., description="该样本泳道中的条带数")


class BandOut(BaseModel):
    center: float = Field(..., description="条带中心长度 (bp)，成员长度的均值")
    min_length: float
    max_length: float
    reference_count: int = Field(..., description="参考型落入该条带的片段数")
    variant_count: int = Field(..., description="变异型落入该条带的片段数")
    reference_lengths: list[float]
    variant_lengths: list[float]
    distinguishable: bool = Field(..., description="两样本条带计数是否不同")


class DifferenceOut(BaseModel):
    band_center: float
    band_range: list[float] = Field(..., description="[min_length, max_length]")
    reference_count: int
    variant_count: int
    reference_lengths: list[float]
    variant_lengths: list[float]
    reason: str = Field(..., description="差异依据说明")


class RecommendationOut(BaseModel):
    rank: int
    enzymes: list[str]
    discriminates: bool = Field(..., description="是否能区分参考型与变异型")
    distinguishable_bands: int = Field(..., description="可区分条带数")
    total_bands: int = Field(..., description="总条带数（两样本合并聚类后的条带数）")
    reference: SampleDigestOut
    variant: SampleDigestOut
    bands: list[BandOut]
    difference_basis: list[DifferenceOut] = Field(..., description="差异依据列表")


class SkippedPairOut(BaseModel):
    enzymes: list[str]
    reason: str


class SummaryOut(BaseModel):
    enzymes: list[str] = Field(..., description="参与枚举的全部候选酶")
    enzyme_pair_count: int = Field(..., description="枚举出的酶对总数")
    evaluated_pairs: int = Field(..., description="通过校验并完成评估的酶对数")
    skipped_pairs: int = Field(..., description="因有效切点不足被拒绝的酶对数")
    discriminating_pairs: int = Field(..., description="能区分两种样本的酶对数")
    returned: int = Field(..., description="本次返回的推荐方案数")


class DigestResponse(BaseModel):
    plasmid_length: int
    min_visible_fragment: float
    band_tolerance: float
    summary: SummaryOut
    skipped: list[SkippedPairOut]
    recommendations: list[RecommendationOut]


class HealthOut(BaseModel):
    status: str
    version: str
