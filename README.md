# 圆形 DNA 双酶切鉴别服务

基于 FastAPI 的环状质粒双酶切方案推荐服务。客户端提交质粒长度、参考型与变异型的
候选酶切位点、最小可见片段与条带分辨容差；服务枚举全部酶对，按圆环相邻切点计算
两种样本的片段长度，将容差内的片段聚为同一条带，并按「可区分条带数 → 总条带数 →
酶名」排序推荐最能区分两种样本的酶切方案。

核心计算为纯 Python 实现，不依赖任何外部服务（无需数据库 / 网络资源）。

## 快速启动（单条命令）

```bash
docker compose up --build
```

后台运行加 `-d`。启动后：

| 地址 | 说明 |
| --- | --- |
| http://localhost:8000 | 服务信息 |
| http://localhost:8000/docs | Swagger 交互式 API 文档 |
| http://localhost:8000/health | 健康检查（`GET`，返回 `{"status":"ok"}`） |
| http://localhost:8000/api/v1/digest/recommend | 酶切方案推荐（`POST`） |

## 配置

均通过环境变量设置，均有默认值，非必配：

| 变量 | 默认值 | 说明 |
| --- | --- | --- |
| `PORT` | `8000` | 宿主机映射端口（容器内固定为 8000），如 `PORT=9000 docker compose up --build` |
| `UVICORN_WORKERS` | `1` | uvicorn 工作进程数 |

## API

### `GET /health`

健康检查，返回 `200 {"status":"ok","version":"1.0.0"}`。

### `POST /api/v1/digest/recommend`

请求体：

| 字段 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `plasmid_length` | int > 0 | 是 | 环状质粒总长度 (bp) |
| `reference_sites` | object | 是 | 参考型候选酶切位点：`{"酶名": [切点...]}` |
| `variant_sites` | object | 是 | 变异型候选酶切位点，结构同上 |
| `min_visible_fragment` | number ≥ 0 | 否（默认 0） | 最小可见片段 (bp)，短于该值的片段不参与条带聚类 |
| `band_tolerance` | number ≥ 0 | 否（默认 0） | 条带分辨容差 (bp)，长度差不超过该值的片段聚为同一条带 |
| `top_n` | int 1–100 | 否（默认全部） | 最多返回的推荐方案数 |

- 切点位置取值范围为 `[0, plasmid_length)`；同一酶对下两种酶的切点取并集。
- 两个样本中未列出的酶按「无切点」处理（可用于模拟变异导致位点丢失）。

示例：

```bash
curl -X POST http://localhost:8000/api/v1/digest/recommend \
  -H 'Content-Type: application/json' \
  -d '{
    "plasmid_length": 6000,
    "reference_sites": {"EcoRI": [500, 3000], "BamHI": [1200], "HindIII": [4000]},
    "variant_sites":  {"EcoRI": [500, 3000], "BamHI": [1200, 4500], "HindIII": [4000]},
    "min_visible_fragment": 0,
    "band_tolerance": 100,
    "top_n": 1
  }'
```

响应（节选）：每个推荐方案包含酶对、两种样本的切点与片段（含 `cross_origin_fragment`
跨原点片段）、聚类条带 `bands`、差异依据 `difference_basis`，以及排序所用的
`distinguishable_bands` / `total_bands`：

```json
{
  "summary": {"enzyme_pair_count": 3, "evaluated_pairs": 3, "skipped_pairs": 0,
              "discriminating_pairs": 2, "returned": 1},
  "recommendations": [
    {
      "rank": 1,
      "enzymes": ["BamHI", "HindIII"],
      "discriminates": true,
      "distinguishable_bands": 3,
      "total_bands": 3,
      "reference": {
        "cut_positions": [1200, 4000],
        "cross_origin_fragment": {"start": 4000, "end": 1200, "length": 3200,
                                  "crosses_origin": true, "visible": true},
        "fragments": ["..."]
      },
      "variant": {"...": "..."},
      "bands": ["..."],
      "difference_basis": [
        {"band_center": 500.0, "reference_count": 0, "variant_count": 1,
         "reason": "约 500.0 bp 处：变异型有 1 条条带，参考型仅 0 条……"}
      ]
    }
  ]
}
```

## 计算规则

- **片段计算**：切点升序后，圆环上相邻切点之差即为片段长度；末切点绕回首切点的
  片段为跨原点片段（`crosses_origin=true`，此时 `end < start`）。
- **条带聚类**：两种样本的可见片段合并排序，相邻长度差 ≤ `band_tolerance` 的片段
  归入同一条带（单链聚类），条带中心取成员均值。
- **可区分条带**：同一条带中参考型与变异型的片段数不同，即为可区分条带，
  全部可区分条带构成 `difference_basis`（差异依据）。
- **推荐排序**：可区分条带数降序 → 总条带数升序（图谱越简单越好）→ 酶名升序。

## 校验与拒绝（422）

- 位点越界：切点不在 `[0, plasmid_length)` 内；
- 重复位点：同一样本同一酶下出现重复切点；
- 候选酶不足两种，无法枚举酶对；
- 任一样本有效切点不足两个：该酶对被拒绝并记录于响应 `skipped`；
  若所有酶对均被拒绝，则整个请求返回 422。

## 本地开发与测试

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements-dev.txt
pytest                      # 运行测试
uvicorn app.main:app --reload --port 8000   # 本地启动
```

## 项目结构

```
app/
├── main.py      # FastAPI 路由：/health、/api/v1/digest/recommend
├── digest.py    # 核心算法：片段计算、条带聚类、酶对枚举与排序（纯 Python）
└── schemas.py   # 请求/响应模型与领域校验
tests/           # pytest 单元与 API 集成测试
Dockerfile
docker-compose.yml
requirements.txt / requirements-dev.txt
```
