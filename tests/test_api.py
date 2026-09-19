"""API 集成测试：推荐排序、跨原点片段、聚类条带、差异依据与拒绝场景。"""
import pytest
from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

# 质粒 6000 bp，三种候选酶；变异型在 BamHI 上多一个切点、HindIII 位点不变。
BASE_PAYLOAD = {
    "plasmid_length": 6000,
    "reference_sites": {
        "EcoRI": [500, 3000],
        "BamHI": [1200],
        "HindIII": [4000],
    },
    "variant_sites": {
        "EcoRI": [500, 3000],
        "BamHI": [1200, 4500],
        "HindIII": [4000],
    },
    "min_visible_fragment": 0,
    "band_tolerance": 100,
}


def post(payload):
    return client.post("/api/v1/digest/recommend", json=payload)


def test_health():
    resp = client.get("/health")
    assert resp.status_code == 200
    assert resp.json()["status"] == "ok"


def test_root():
    resp = client.get("/")
    assert resp.status_code == 200
    assert resp.json()["recommend"] == "/api/v1/digest/recommend"


def test_recommendation_ranking():
    resp = post(BASE_PAYLOAD)
    assert resp.status_code == 200
    body = resp.json()

    assert body["summary"]["enzyme_pair_count"] == 3
    assert body["summary"]["evaluated_pairs"] == 3
    assert body["summary"]["skipped_pairs"] == 0
    assert body["summary"]["discriminating_pairs"] == 2

    recs = body["recommendations"]
    assert [r["rank"] for r in recs] == [1, 2, 3]

    # 可区分条带数降序 → 总条带数升序 → 酶名升序
    first, second, third = recs
    assert first["enzymes"] == ["BamHI", "HindIII"]
    assert first["distinguishable_bands"] == 3
    assert first["total_bands"] == 3
    assert second["enzymes"] == ["BamHI", "EcoRI"]
    assert second["distinguishable_bands"] == 3
    assert second["total_bands"] == 5
    assert third["enzymes"] == ["EcoRI", "HindIII"]
    assert third["distinguishable_bands"] == 0
    assert third["discriminates"] is False
    assert third["difference_basis"] == []


def test_fragments_and_cross_origin():
    resp = post(BASE_PAYLOAD)
    rec = next(
        r for r in resp.json()["recommendations"] if r["enzymes"] == ["BamHI", "EcoRI"]
    )

    ref = rec["reference"]
    assert ref["cut_positions"] == [500, 1200, 3000]
    assert sorted(f["length"] for f in ref["fragments"]) == [700, 1800, 3500]
    cross = ref["cross_origin_fragment"]
    assert cross["crosses_origin"] is True
    assert (cross["start"], cross["end"], cross["length"]) == (3000, 500, 3500)

    var = rec["variant"]
    assert var["cut_positions"] == [500, 1200, 3000, 4500]
    assert sorted(f["length"] for f in var["fragments"]) == [700, 1500, 1800, 2000]


def test_band_clustering_and_difference_basis():
    resp = post(BASE_PAYLOAD)
    rec = next(
        r for r in resp.json()["recommendations"] if r["enzymes"] == ["BamHI", "EcoRI"]
    )

    bands = rec["bands"]
    assert len(bands) == 5
    by_center = {round(b["center"]): b for b in bands}
    # 共有条带：700 / 1800，两样本各一条，不可区分
    assert by_center[700]["reference_count"] == 1
    assert by_center[700]["variant_count"] == 1
    assert by_center[700]["distinguishable"] is False
    # 差异条带：1500、2000 仅变异型；3500 仅参考型
    assert by_center[1500]["variant_count"] == 1
    assert by_center[1500]["reference_count"] == 0
    assert by_center[3500]["reference_count"] == 1
    assert by_center[3500]["variant_count"] == 0

    basis = rec["difference_basis"]
    assert len(basis) == 3
    assert {round(d["band_center"]) for d in basis} == {1500, 2000, 3500}
    assert all(d["reason"] for d in basis)


def test_tolerance_merges_close_bands():
    resp = post(BASE_PAYLOAD)
    rec = next(
        r
        for r in resp.json()["recommendations"]
        if r["enzymes"] == ["BamHI", "HindIII"]
    )
    # 变异型 2700 与共有 2800 相差 100，容差 100 内聚为同一条带
    merged = next(b for b in rec["bands"] if b["min_length"] == 2700)
    assert merged["max_length"] == 2800
    assert merged["reference_count"] == 1
    assert merged["variant_count"] == 2
    assert merged["distinguishable"] is True


def test_min_visible_fragment_filters_bands():
    payload = dict(BASE_PAYLOAD, min_visible_fragment=600)
    resp = post(payload)
    assert resp.status_code == 200
    rec = next(
        r
        for r in resp.json()["recommendations"]
        if r["enzymes"] == ["BamHI", "HindIII"]
    )
    var = rec["variant"]
    invisible = [f for f in var["fragments"] if not f["visible"]]
    assert [f["length"] for f in invisible] == [500]
    assert 500 not in var["visible_fragment_lengths"]
    # 500 bp 片段不可见后，该酶对只剩 2 个可区分条带
    assert rec["distinguishable_bands"] == 2


def test_top_n_limits_recommendations():
    resp = post(dict(BASE_PAYLOAD, top_n=1))
    assert resp.status_code == 200
    body = resp.json()
    assert len(body["recommendations"]) == 1
    assert body["summary"]["returned"] == 1
    assert body["summary"]["evaluated_pairs"] == 3


def test_out_of_bounds_site_rejected():
    payload = dict(BASE_PAYLOAD)
    payload["reference_sites"] = dict(BASE_PAYLOAD["reference_sites"], EcoRI=[500, 6000])
    resp = post(payload)
    assert resp.status_code == 422
    assert "越界" in resp.text


def test_negative_site_rejected():
    payload = dict(BASE_PAYLOAD)
    payload["variant_sites"] = dict(BASE_PAYLOAD["variant_sites"], HindIII=[-1])
    resp = post(payload)
    assert resp.status_code == 422


def test_duplicate_site_rejected():
    payload = dict(BASE_PAYLOAD)
    payload["reference_sites"] = dict(BASE_PAYLOAD["reference_sites"], EcoRI=[500, 500])
    resp = post(payload)
    assert resp.status_code == 422
    assert "重复位点" in resp.text


def test_single_enzyme_rejected():
    resp = post(
        {
            "plasmid_length": 6000,
            "reference_sites": {"EcoRI": [500, 3000]},
            "variant_sites": {"EcoRI": [500, 3000]},
        }
    )
    assert resp.status_code == 422


def test_insufficient_cuts_all_pairs_rejected():
    resp = post(
        {
            "plasmid_length": 6000,
            "reference_sites": {"EcoRI": [500], "BamHI": []},
            "variant_sites": {"EcoRI": [500], "BamHI": []},
        }
    )
    assert resp.status_code == 422
    assert "有效切点不足" in resp.text


def test_pair_with_insufficient_cuts_skipped():
    resp = post(
        {
            "plasmid_length": 6000,
            "reference_sites": {
                "EcoRI": [500, 3000],
                "BamHI": [1200],
                "HindIII": [4000],
            },
            "variant_sites": {
                "EcoRI": [500, 3000],
                "BamHI": [1200],
                "HindIII": [],
            },
        }
    )
    assert resp.status_code == 200
    body = resp.json()
    # (BamHI, HindIII)：变异型合并切点仅 1 个 → 该酶对被拒绝
    assert body["summary"]["skipped_pairs"] == 1
    assert body["summary"]["evaluated_pairs"] == 2
    assert body["skipped"][0]["enzymes"] == ["BamHI", "HindIII"]
    assert "有效切点不足" in body["skipped"][0]["reason"]


def test_invalid_plasmid_length_rejected():
    resp = post(dict(BASE_PAYLOAD, plasmid_length=0))
    assert resp.status_code == 422


def test_negative_tolerance_rejected():
    resp = post(dict(BASE_PAYLOAD, band_tolerance=-5))
    assert resp.status_code == 422
