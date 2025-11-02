# GIM3E: Gene Inactivation Moderated by Metabolism, Metabolomics, and Expression

[English](#english) | [한국어](#korean)

---

<a name="english"></a>
## English Documentation

### Overview

**GIM3E** (Gene Inactivity Moderated by Metabolism, Metabolomics, and Expression) is a computational framework for integrating multi-omics data with genome-scale metabolic models to create condition-specific metabolic models.

**Key Capabilities:**
- 🧬 **Transcriptomics Integration**: Incorporates gene expression data (RNA-seq, microarray, proteomics)
- 🔬 **Metabolomics Integration**: Uses detected metabolites to constrain metabolic flux
- 🧮 **Constraint-based Modeling**: Built on COBRApy framework
- 📊 **Flux Variability Analysis**: Determines flux ranges for all reactions
- 🎯 **Condition-specific Models**: Creates tailored models for specific biological conditions

### Citation

If you use GIM3E in your research, please cite:

> Schmidt BJ, Ebrahim A, Metz TO, Adkins JN, Palsson BØ, Hyduke DR. (2013)
> **GIM3E: condition-specific models of cellular metabolism developed from metabolomics and expression data.**
> *Bioinformatics*, 29(22):2900-2908.
> doi: [10.1093/bioinformatics/btt493](https://doi.org/10.1093/bioinformatics/btt493)

### Installation

#### Prerequisites

**Python Version:** Python 3.7 or higher

**Required Dependencies:**
- COBRApy >= 0.26.0
- NumPy >= 1.20.0
- SciPy >= 1.7.0
- Pandas >= 1.3.0

**Recommended Solvers:**
- **CPLEX** 12.10+ (Commercial, IBM Academic License available) - **Strongly Recommended**
- **Gurobi** 9.0+ (Commercial, Academic License available)
- GLPK (Open-source, included with COBRApy, but poor MILP performance)

#### Installation Steps

1. **Clone the repository:**
   ```bash
   git clone https://github.com/jyryu3161/gim3e.git
   cd gim3e
   ```

2. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

3. **Install GIM3E:**
   ```bash
   pip install -e .
   ```

4. **Install a solver (choose one):**

   **Option A: CPLEX (Recommended)**
   ```bash
   pip install cplex
   ```

   **Option B: Gurobi**
   ```bash
   pip install gurobipy
   ```

   **Option C: GLPK (Default)**
   - Automatically installed with COBRApy

### Quick Start

#### Basic Example

```python
from gim3e.core import gim3e
from cobra.io import load_model

# 1. Load a metabolic model
model = load_model("textbook")  # E. coli core model

# 2. Prepare transcriptomics data
expression_dict = {
    'b0116': 8.5,  # High expression
    'b0118': 3.2,  # Low expression
    'b0734': 6.7,
    'b0733': 7.1,
}

# 3. Prepare metabolomics data
detected_metabolites = [
    'glc__D_c',   # Glucose
    'pyr_c',      # Pyruvate
    'lac__D_c',   # Lactate
    'atp_c',      # ATP
]

# 4. Run GIM3E
gim3e_model, fva_results, penalty_score = gim3e.gim3e(
    model,
    expression_dict=expression_dict,
    expression_threshold=5.0,
    metabolite_list=detected_metabolites,
    fraction_growth=0.9,
    MILP_formulation=True,
    run_FVA=True,
    solver='cplex'
)

# 5. Analyze results
print(f"Total penalty score: {penalty_score}")
print(f"Active reactions: {len([r for r in fva_results if fva_results[r]['maximum'] > 1e-6])}")
```

#### Recon 2 Human Model Example

See the comprehensive example:
```bash
cd gim3e/examples
python 03_recon2_metabolomics_transcriptomics_integration.py
```

This example demonstrates:
- Loading Recon 2 human metabolic model
- Simulating cancer vs normal cell metabolism
- Integrating RNA-seq and LC-MS data
- Analyzing condition-specific metabolic changes

### Examples

The `gim3e/examples/` directory contains:

1. **01_run_gim3e_with_salmonella.py**
   - Original Salmonella enterica example
   - Full GIM3E workflow with validation

2. **02_sample_core_gim3e.py**
   - E. coli core model with GMS sampling

3. **03_recon2_metabolomics_transcriptomics_integration.py** ⭐ **NEW**
   - Human Recon 2 model example
   - Cancer vs normal cell metabolism
   - Modern COBRApy API usage
   - Complete workflow from data to results

### Key Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cobra_model` | Model | Required | COBRApy metabolic model |
| `expression_dict` | dict | `{}` | Gene expression data {gene_id: value} |
| `expression_threshold` | float | 0.0 | Cutoff for applying expression penalties |
| `metabolite_list` | list | `[]` | Detected metabolite IDs |
| `fraction_growth` | float | 0.9 | Minimum fraction of optimal growth (0-1) |
| `relative_penalty_bound` | float | 1.0 | Maximum penalty increase allowed (≥1.0) |
| `MILP_formulation` | bool | `False` | Use MILP to eliminate flux loops |
| `run_FVA` | bool | `True` | Perform flux variability analysis |
| `solver` | str | `'cplex'` | Solver to use (cplex, gurobi, glpk) |

### Troubleshooting

**Solver not found:**
```bash
pip install cplex  # or gurobipy
# or use solver='glpk' (not recommended for MILP)
```

**Infeasible model:**
- Check metabolite IDs match model exactly
- Verify model optimizes: `model.optimize()`
- Reduce `fraction_growth` (e.g., from 0.9 to 0.7)

**Slow optimization:**
- Use CPLEX or Gurobi instead of GLPK
- Set `MILP_formulation=False` for faster LP
- Reduce model size or metabolite_list

### Performance

| Model Size | Reactions | MILP Time (CPLEX) | LP Time (CPLEX) |
|------------|-----------|-------------------|-----------------|
| E. coli core | ~100 | < 1 min | < 10 sec |
| E. coli iJO1366 | ~2,500 | 10-30 min | 1-5 min |
| Recon 2 | ~7,500 | 1-3 hours | 5-15 min |
| Recon 3D | ~10,000+ | 2-6 hours | 10-30 min |

### License

GIM3E is licensed under the **GNU General Public License v3.0**.

### Support

- **Issues:** https://github.com/jyryu3161/gim3e/issues
- **Paper:** https://doi.org/10.1093/bioinformatics/btt493
- **COBRApy Docs:** https://cobrapy.readthedocs.io/

---

<a name="korean"></a>
## 한국어 문서

### 개요

**GIM3E** (Gene Inactivity Moderated by Metabolism, Metabolomics, and Expression)는 다중 오믹스 데이터를 전장 유전체 규모 대사 모델과 통합하여 조건 특이적 대사 모델을 생성하는 계산 프레임워크입니다.

**주요 기능:**
- 🧬 **전사체 통합**: 유전자 발현 데이터 통합 (RNA-seq, 마이크로어레이, 프로테오믹스)
- 🔬 **대사체 통합**: 검출된 대사산물을 사용하여 대사 플럭스 제약
- 🧮 **제약 기반 모델링**: COBRApy 프레임워크 기반
- 📊 **플럭스 가변성 분석**: 모든 반응의 플럭스 범위 결정
- 🎯 **조건 특이적 모델**: 특정 생물학적 조건에 맞춤화된 모델 생성

### 인용

연구에서 GIM3E를 사용하는 경우 다음을 인용해 주세요:

> Schmidt BJ, Ebrahim A, Metz TO, Adkins JN, Palsson BØ, Hyduke DR. (2013)
> **GIM3E: condition-specific models of cellular metabolism developed from metabolomics and expression data.**
> *Bioinformatics*, 29(22):2900-2908.
> doi: [10.1093/bioinformatics/btt493](https://doi.org/10.1093/bioinformatics/btt493)

### 설치

#### 필수 요구사항

**Python 버전:** Python 3.7 이상

**필수 의존성:**
- COBRApy >= 0.26.0
- NumPy >= 1.20.0
- SciPy >= 1.7.0
- Pandas >= 1.3.0

**권장 솔버:**
- **CPLEX** 12.10+ (상용, IBM 학술 라이선스 이용 가능) - **강력 권장**
- **Gurobi** 9.0+ (상용, 학술 라이선스 이용 가능)
- GLPK (오픈소스, COBRApy에 포함, MILP 성능 낮음)

#### 설치 단계

1. **저장소 복제:**
   ```bash
   git clone https://github.com/jyryu3161/gim3e.git
   cd gim3e
   ```

2. **의존성 설치:**
   ```bash
   pip install -r requirements.txt
   ```

3. **GIM3E 설치:**
   ```bash
   pip install -e .
   ```

4. **솔버 설치 (하나 선택):**

   **옵션 A: CPLEX (권장)**
   ```bash
   pip install cplex
   ```

   **옵션 B: Gurobi**
   ```bash
   pip install gurobipy
   ```

   **옵션 C: GLPK (기본)**
   - COBRApy와 함께 자동 설치

### 빠른 시작

#### 기본 예제

```python
from gim3e.core import gim3e
from cobra.io import load_model

# 1. 대사 모델 로드
model = load_model("textbook")  # E. coli 핵심 모델

# 2. 전사체 데이터 준비
expression_dict = {
    'b0116': 8.5,  # 높은 발현
    'b0118': 3.2,  # 낮은 발현
    'b0734': 6.7,
    'b0733': 7.1,
}

# 3. 대사체 데이터 준비
detected_metabolites = [
    'glc__D_c',   # 포도당
    'pyr_c',      # 피루브산
    'lac__D_c',   # 젖산
    'atp_c',      # ATP
]

# 4. GIM3E 실행
gim3e_model, fva_results, penalty_score = gim3e.gim3e(
    model,
    expression_dict=expression_dict,
    expression_threshold=5.0,
    metabolite_list=detected_metabolites,
    fraction_growth=0.9,
    MILP_formulation=True,
    run_FVA=True,
    solver='cplex'
)

# 5. 결과 분석
print(f"총 페널티 점수: {penalty_score}")
print(f"활성 반응 수: {len([r for r in fva_results if fva_results[r]['maximum'] > 1e-6])}")
```

#### Recon 2 인간 모델 예제

전체 예제 참조:
```bash
cd gim3e/examples
python 03_recon2_metabolomics_transcriptomics_integration.py
```

이 예제는 다음을 보여줍니다:
- Recon 2 인간 대사 모델 로딩
- 암세포 vs 정상세포 대사 시뮬레이션
- RNA-seq 및 LC-MS 데이터 통합
- 조건 특이적 대사 변화 분석

### 예제

`gim3e/examples/` 디렉토리에는 다음이 포함됩니다:

1. **01_run_gim3e_with_salmonella.py**
   - 원본 Salmonella enterica 예제
   - 검증이 포함된 전체 GIM3E 워크플로우

2. **02_sample_core_gim3e.py**
   - GMS 샘플링이 포함된 E. coli 핵심 모델

3. **03_recon2_metabolomics_transcriptomics_integration.py** ⭐ **신규**
   - 인간 Recon 2 모델 예제
   - 암세포 vs 정상세포 대사
   - 최신 COBRApy API 사용
   - 데이터에서 결과까지 완전한 워크플로우

### 주요 매개변수

| 매개변수 | 타입 | 기본값 | 설명 |
|---------|------|--------|------|
| `cobra_model` | Model | 필수 | COBRApy 대사 모델 |
| `expression_dict` | dict | `{}` | 유전자 발현 데이터 {gene_id: value} |
| `expression_threshold` | float | 0.0 | 발현 페널티 적용 임계값 |
| `metabolite_list` | list | `[]` | 검출된 대사산물 ID |
| `fraction_growth` | float | 0.9 | 최적 성장의 최소 비율 (0-1) |
| `relative_penalty_bound` | float | 1.0 | 허용되는 최대 페널티 증가 (≥1.0) |
| `MILP_formulation` | bool | `False` | MILP를 사용하여 플럭스 루프 제거 |
| `run_FVA` | bool | `True` | 플럭스 가변성 분석 수행 |
| `solver` | str | `'cplex'` | 사용할 솔버 (cplex, gurobi, glpk) |

### 문제 해결

**솔버를 찾을 수 없음:**
```bash
pip install cplex  # 또는 gurobipy
# 또는 solver='glpk' 사용 (MILP에는 권장하지 않음)
```

**실행 불가능한 모델:**
- 대사산물 ID가 모델과 정확히 일치하는지 확인
- 모델 최적화 확인: `model.optimize()`
- `fraction_growth` 감소 (예: 0.9에서 0.7로)

**느린 최적화:**
- GLPK 대신 CPLEX 또는 Gurobi 사용
- 더 빠른 LP를 위해 `MILP_formulation=False` 설정
- 모델 크기 또는 metabolite_list 감소

### 성능

| 모델 크기 | 반응 수 | MILP 시간 (CPLEX) | LP 시간 (CPLEX) |
|----------|--------|------------------|----------------|
| E. coli core | ~100 | < 1분 | < 10초 |
| E. coli iJO1366 | ~2,500 | 10-30분 | 1-5분 |
| Recon 2 | ~7,500 | 1-3시간 | 5-15분 |
| Recon 3D | ~10,000+ | 2-6시간 | 10-30분 |

### 라이선스

GIM3E는 **GNU General Public License v3.0**으로 라이선스가 부여됩니다.

### 지원

- **이슈:** https://github.com/jyryu3161/gim3e/issues
- **논문:** https://doi.org/10.1093/bioinformatics/btt493
- **COBRApy 문서:** https://cobrapy.readthedocs.io/

---

## Version History

### v2.0.0 (2025)
- ✨ Updated for modern Python 3.7+ and COBRApy 0.26+
- 📝 Added comprehensive bilingual README (English/Korean)
- 🔬 Added Recon 2 human metabolic model integration example
- 🐛 Modernized setup.py and dependencies
- 📦 Added requirements.txt
- 🎯 Added detailed documentation

### v1.0.3 (2013)
- Substantial GMS improvements
- Original release with Python 2.7 support

---

**Original Development:** University of California, San Diego, Systems Biology Research Group
**Funding:** NIAID interagency agreement Y1-AI-8401
**Maintainer:** Updated for modern Python ecosystem
**Original Author:** Brian J Schmidt
**Last Updated:** 2025-11-02
