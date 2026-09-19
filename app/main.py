"""FastAPI 入口：圆形 DNA 双酶切鉴别服务。"""
from __future__ import annotations

from fastapi import FastAPI, HTTPException

from app import __version__
from app.digest import Band, PairResult, RecommendationResult, recommend
from app.schemas import (
    BandOut,
    DifferenceOut,
    DigestRequest,
    DigestResponse,
    FragmentOut,
    HealthOut,
    RecommendationOut,
    SampleDigestOut,
    SkippedPairOut,
    SummaryOut,
)

app = FastAPI(
    title="圆形 DNA 双酶切鉴别服务",
    description=(
        "提交质粒长度、参考型/变异型候选酶切位点、最小可见片段与条带分辨容差，"
        "服务枚举酶对、计算环状片段、聚类条带并推荐可区分两种样本的酶切方案。"
    ),
    version=__version__,
)


@app.get("/health", response_model=HealthOut, tags=["meta"], summary="健康检查")
def health() -> HealthOut:
    return HealthOut(status="ok", version=__version__)


@app.get("/", tags=["meta"], summary="服务信息")
def root() -> dict:
    return {
        "service": "circular-dna-double-digest",
        "version": __version__,
        "docs": "/docs",
        "health": "/health",
        "recommend": "/api/v1/digest/recommend",
    }


def _fragment_out(fragment) -> FragmentOut:
    return FragmentOut(
        start=fragment.start,
        end=fragment.end,
        length=fragment.length,
        crosses_origin=fragment.crosses_origin,
        visible=fragment.visible,
    )


def _sample_out(sample) -> SampleDigestOut:
    cross_origin = next(f for f in sample.fragments if f.crosses_origin)
    return SampleDigestOut(
        cut_positions=sample.cut_positions,
        fragment_count=len(sample.fragments),
        fragments=[_fragment_out(f) for f in sample.fragments],
        cross_origin_fragment=_fragment_out(cross_origin),
        visible_fragment_lengths=[f.length for f in sample.fragments if f.visible],
        band_count=sample.band_count,
    )


def _band_out(band: Band) -> BandOut:
    return BandOut(
        center=round(band.center, 3),
        min_length=band.min_length,
        max_length=band.max_length,
        reference_count=band.reference_count,
        variant_count=band.variant_count,
        reference_lengths=band.reference_lengths,
        variant_lengths=band.variant_lengths,
        distinguishable=band.distinguishable,
    )


def _difference_out(band: Band) -> DifferenceOut:
    center = round(band.center, 3)
    if band.reference_count > band.variant_count:
        reason = (
            f"约 {center} bp 处：参考型有 {band.reference_count} 条条带，"
            f"变异型仅 {band.variant_count} 条（变异型切点变化导致该片段丢失或合并）"
        )
    else:
        reason = (
            f"约 {center} bp 处：变异型有 {band.variant_count} 条条带，"
            f"参考型仅 {band.reference_count} 条（变异型新增/丢失切点导致片段数目变化）"
        )
    return DifferenceOut(
        band_center=center,
        band_range=[band.min_length, band.max_length],
        reference_count=band.reference_count,
        variant_count=band.variant_count,
        reference_lengths=band.reference_lengths,
        variant_lengths=band.variant_lengths,
        reason=reason,
    )


def _recommendation_out(rank: int, pair: PairResult) -> RecommendationOut:
    return RecommendationOut(
        rank=rank,
        enzymes=list(pair.enzymes),
        discriminates=pair.discriminates,
        distinguishable_bands=pair.distinguishable_bands,
        total_bands=pair.total_bands,
        reference=_sample_out(pair.reference),
        variant=_sample_out(pair.variant),
        bands=[_band_out(b) for b in pair.bands],
        difference_basis=[_difference_out(b) for b in pair.bands if b.distinguishable],
    )


@app.post(
    "/api/v1/digest/recommend",
    response_model=DigestResponse,
    tags=["digest"],
    summary="枚举酶对并推荐可区分参考型/变异型的双酶切方案",
)
def recommend_endpoint(request: DigestRequest) -> DigestResponse:
    result: RecommendationResult = recommend(
        plasmid_length=request.plasmid_length,
        reference_sites=request.reference_sites,
        variant_sites=request.variant_sites,
        min_visible_fragment=request.min_visible_fragment,
        band_tolerance=request.band_tolerance,
    )
    if not result.pairs:
        raise HTTPException(
            status_code=422,
            detail=(
                "任一样本有效切点不足两个，所有酶对均被拒绝，"
                "无法进行双酶切分析"
            ),
        )

    pairs = result.pairs[: request.top_n] if request.top_n else result.pairs
    return DigestResponse(
        plasmid_length=request.plasmid_length,
        min_visible_fragment=request.min_visible_fragment,
        band_tolerance=request.band_tolerance,
        summary=SummaryOut(
            enzymes=sorted(set(request.reference_sites) | set(request.variant_sites)),
            enzyme_pair_count=result.enzyme_pair_count,
            evaluated_pairs=len(result.pairs),
            skipped_pairs=len(result.skipped),
            discriminating_pairs=sum(1 for p in result.pairs if p.discriminates),
            returned=len(pairs),
        ),
        skipped=[
            SkippedPairOut(enzymes=list(s.enzymes), reason=s.reason)
            for s in result.skipped
        ],
        recommendations=[
            _recommendation_out(rank, pair) for rank, pair in enumerate(pairs, start=1)
        ],
    )
