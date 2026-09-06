# mapping.py
"""
MAPS III 2025 胃镜监测决策模块（符合新 Schema）
"""
from collections import OrderedDict
import re
from enum import Enum
from typing import Dict, Any, List, Optional, Union

# =============================================================================
# DECISION STATUS
# =============================================================================
class DecisionStatus(str, Enum):
    TRUE = "TRUE"
    FALSE = "FALSE"
    UNKNOWN = "UNKNOWN"


# =============================================================================
# RESULT FACTORY
# =============================================================================
def make_result(
    status: DecisionStatus,
    missing_fields: Optional[List[str]] = None,
    extra: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    result: Dict[str, Any] = {"status": status}
    if missing_fields:
        result["missing_fields"] = sorted(list(set(missing_fields)))
    if extra:
        result.update(extra)
    return result


# =============================================================================
# BASIC HELPERS
# =============================================================================
def require_fields(data: Dict[str, Any], fields: List[str]) -> List[str]:
    """返回缺失（不存在或为 None 或 'not_mentioned'）的字段列表"""
    missing = []
    for f in fields:
        val = data.get(f)
        if val is None or val == "not_mentioned":
            missing.append(f)
    return missing

def _safe_float(x: Any) -> Optional[float]:
    try:
        if x is None:
            return None
        if isinstance(x, str):
            match = re.search(r'(\d+(?:\.\d+)?)', x)
            if match:
                return float(match.group(1))
            else:
                return None
        return float(x)
    except (TypeError, ValueError):
        return None

def _safe_int(x: Any) -> Optional[int]:
    try:
        if x is None:
            return None
        if isinstance(x, str):
            match = re.search(r'(\d+)', x)
            if match:
                return int(match.group(1))
            else:
                return None
        return int(x)
    except (TypeError, ValueError):
        return None

def _parse_measurement(value: Any) -> Optional[float]:
    return _safe_float(value)

def _yes_no_to_bool(val: Any) -> Optional[bool]:
    if val is None or val == "not_mentioned":
        return None
    if isinstance(val, str):
        if val.lower() == "yes":
            return True
        elif val.lower() == "no":
            return False
    return None

def _enum_match(val: Any, expected: str) -> bool:
    return val == expected


# =============================================================================
# PEPSINOGEN
# =============================================================================
def pepsinogen_ratio(data: Dict[str, Any]) -> Optional[float]:
    pgi = _safe_float(data.get("pepsinogen_i_level"))
    pgii = _safe_float(data.get("pepsinogen_ii_level"))
    if pgi is None or pgii is None or pgii == 0:
        return None
    return pgi / pgii

def high_serologic_risk(data: Dict[str, Any]) -> Optional[bool]:
    pgi = _safe_float(data.get("pepsinogen_i_level"))
    ratio = pepsinogen_ratio(data)
    if pgi is None and ratio is None:
        return None
    cond_pgi = pgi is not None and pgi <= 70
    cond_ratio = ratio is not None and ratio <= 3
    return cond_pgi or cond_ratio


# =============================================================================
# 独立特征判定函数（精确支持“或”逻辑）
# =============================================================================
def _kt_high(data: Dict[str, Any]) -> Optional[bool]:
    val = data.get("kimura_takemoto_classification")
    if val is None or val == "not_mentioned":
        return None
    return val.lower() in ["c3", "o1", "o2", "o3"]

def _kt_low(data: Dict[str, Any]) -> Optional[bool]:
    val = data.get("kimura_takemoto_classification")
    if val is None or val == "not_mentioned":
        return None
    return val.lower() in ["c0", "c1", "c2"]

def _eggim_high(data: Dict[str, Any]) -> Optional[bool]:
    val = data.get("eggim_score")
    if val is None or val == "not_mentioned":
        return None
    try:
        return int(val) >= 5
    except:
        return None

def _eggim_low(data: Dict[str, Any]) -> Optional[bool]:
    val = data.get("eggim_score")
    if val is None or val == "not_mentioned":
        return None
    try:
        return int(val) <= 4
    except:
        return None

def _olga_high(data: Dict[str, Any]) -> Optional[bool]:
    val = data.get("olga_stage")
    if val is None or val == "not_mentioned":
        return None
    if isinstance(val, str):
        return val.lower() in ["iii", "iv"]
    try:
        return int(val) >= 3
    except:
        return None

def _olga_low(data: Dict[str, Any]) -> Optional[bool]:
    val = data.get("olga_stage")
    if val is None or val == "not_mentioned":
        return None
    if isinstance(val, str):
        return val.lower() in ["0", "i", "ii"]
    try:
        return int(val) <= 2
    except:
        return None

def _olgim_high(data: Dict[str, Any]) -> Optional[bool]:
    val = data.get("olgim_stage")
    if val is None or val == "not_mentioned":
        return None
    if isinstance(val, str):
        return val.lower() in ["iii", "iv"]
    try:
        return int(val) >= 3
    except:
        return None

def _olgim_low(data: Dict[str, Any]) -> Optional[bool]:
    val = data.get("olgim_stage")
    if val is None or val == "not_mentioned":
        return None
    if isinstance(val, str):
        return val.lower() in ["0", "i", "ii"]
    try:
        return int(val) <= 2
    except:
        return None

def has_additional_risk_factor(data: Dict[str, Any]) -> bool:
    if _yes_no_to_bool(data.get("incomplete_intestinal_metaplasia")) is True:
        return True
    if _yes_no_to_bool(data.get("family_history")) is True:
        return True
    hp = data.get("helicobacter_pylori_status")
    if hp == "persistent": 
        return True
    return False


# 旧函数兼容
def high_kimura(data): return _kt_high(data)
def low_kimura(data): return _kt_low(data)
def high_eggim(data): return _eggim_high(data)
def low_eggim(data): return _eggim_low(data)


# =============================================================================
# SM INVASION & EMR ELIGIBLE
# =============================================================================
def sm1_invasion(data: Dict[str, Any]) -> Optional[bool]:
    missing = require_fields(data, ["submucosal_invasion", "submucosal_invasion_depth"])
    if missing:
        return None
    if data["submucosal_invasion"] != "yes":
        return False
    depth = data.get("submucosal_invasion_depth")
    if depth == "le_500um":
        return True
    elif depth == "gt_500um":
        return False
    else:
        return None

def emr_eligible(data: Dict[str, Any]) -> Optional[bool]:
    paris = data.get("lesion_paris_classification")
    if paris != "0_iia":
        return False
    size = _parse_measurement(data.get("neoplastic_lesion_size"))
    if size is None or size > 10:
        return False
    if (_yes_no_to_bool(data.get("neoplastic_lesion_ulceration")) is True or
        data.get("differentiation_status") != "differentiated" or
        data.get("submucosal_invasion") == "yes"):
        return False
    return True
    
# =============================================================================
# NON VISIBLE
# =============================================================================
def _has_nonvisible_dysplasia_or_indefinite(data: Dict[str, Any]) -> Optional[bool]:
    """
    Returns True if there is non‑visible dysplasia (any grade) OR
    indefinite for dysplasia in the absence of a visible lesion.
    """
    dysplasia = data.get("dysplasia")
    indefinite = data.get("indefinite_for_dysplasia")
    grade = data.get("dysplasia_grade")
    visibility = data.get("dysplasia_visibility")
    lesion_vis = data.get("neoplastic_lesion_visibility")

    # If a visible lesion is present, this is NOT a "non‑visible" scenario
    if visibility == "visible" or lesion_vis == "visible":
        return False

    # Positive finding in absence of visible lesion
    positive = (
        dysplasia in ("yes", "highly_likely") or
        _yes_no_to_bool(indefinite) is True or
        grade in ("low_grade", "high_grade")
    )
    required = [dysplasia, indefinite, grade, visibility, lesion_vis]
    if any(v is None or v == "not_mentioned" for v in required):
        return None
    return positive
    
# =============================================================================
# CHAIN 1 – TRULY LOW RISK（需所有低风险字段明确存在）
# =============================================================================
def evaluate_chain_1(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Truly low risk: no extensive atrophy/IM, no risk factors.
    Mild atrophy (K-T C0-2, OLGA 0-II) is acceptable.
    """
    # 必须排除的高危客观指标
    if _kt_high(data) or _eggim_high(data) or _olga_high(data) or _olgim_high(data):
        return make_result(DecisionStatus.FALSE)
    if _yes_no_to_bool(data.get("incomplete_intestinal_metaplasia")) is True:
        return make_result(DecisionStatus.FALSE)

    # 允许轻度萎缩（低 K-T/OLGA）
    # 已由上面排除高危，所以自动允许低危或未提及的萎缩

    # 风险因素（必须明确为 no 或缺失）
    incomplete = data.get("incomplete_intestinal_metaplasia")
    if incomplete not in [None, "not_mentioned"] and incomplete != "no":
        return make_result(DecisionStatus.FALSE)

    family = data.get("family_history")
    if family not in [None, "not_mentioned"] and family != "no":
        return make_result(DecisionStatus.FALSE)

    # 萎缩不单独排除
    # atrophic_gastritis 即使为 yes，只要 K-T/OLGA 低风险，仍可进入 Chain 1
    return make_result(DecisionStatus.TRUE, extra={"chain": 1})


# =============================================================================
# CHAIN 2 – EXTENSIVE IM / HIGH RISK (NO ADDITIONAL RISK)
# =============================================================================
def evaluate_chain_2(data: Dict[str, Any]) -> Dict[str, Any]:
    # 获取所有高危特征
    kt = _kt_high(data)
    eggim = _eggim_high(data)
    olga = _olga_high(data)
    olgim = _olgim_high(data)
    incomplete = _yes_no_to_bool(data.get("incomplete_intestinal_metaplasia"))
    features = [kt, eggim, olga, olgim, incomplete]

    if any(f is True for f in features):
        # 检查是否由 OLGA/OLGIM III/IV 触发 (+REC.37)
        advanced_olg = (_olga_high(data) is True) or (_olgim_high(data) is True)
        return make_result(DecisionStatus.TRUE, extra={"chain": 2, "advanced_olg_staging": advanced_olg})
    if all(f is None for f in features):
        # 完全没有任何风险信息
        return make_result(DecisionStatus.UNKNOWN,
                           ["kimura_takemoto_classification", "eggim_score",
                            "olga_stage", "olgim_stage", "incomplete_intestinal_metaplasia"])
    # 有明确的 False 且无 True → 不是高危
    return make_result(DecisionStatus.FALSE)


# =============================================================================
# CHAIN 3 – EXTENSIVE IM / HIGH RISK + FIRST-DEGREE FAMILY HISTORY
# =============================================================================
def evaluate_chain_3(data: Dict[str, Any]) -> Dict[str, Any]:
    # 先检查高危特征（同Chain 2）
    kt = _kt_high(data)
    eggim = _eggim_high(data)
    olga = _olga_high(data)
    olgim = _olgim_high(data)
    incomplete = _yes_no_to_bool(data.get("incomplete_intestinal_metaplasia"))
    features = [kt, eggim, olga, olgim, incomplete]

    if any(f is True for f in features):
        # 高危存在 → 检查家族史
        family = _yes_no_to_bool(data.get("family_history"))
        if family is True:
            # 检查是否由 OLGA/OLGIM III/IV 触发 (+REC.37)
            advanced_olg = (_olga_high(data) is True) or (_olgim_high(data) is True)
            return make_result(DecisionStatus.TRUE, extra={"chain": 3, "advanced_olg_staging": advanced_olg})
        elif family is False:
            return make_result(DecisionStatus.FALSE)
        else:
            return make_result(DecisionStatus.UNKNOWN, ["family_history"])
    if all(f is None for f in features):
        return make_result(DecisionStatus.UNKNOWN,
                           ["kimura_takemoto_classification", "eggim_score",
                            "olga_stage", "olgim_stage", "incomplete_intestinal_metaplasia"])
    return make_result(DecisionStatus.FALSE)


# =============================================================================
# CHAIN 4 – LIMITED IM + ADDITIONAL RISK FACTOR (per Guideline 35 & Fig.10)
# =============================================================================

def evaluate_chain_4(data: Dict[str, Any]) -> Dict[str, Any]:
    olgim = data.get("olgim_stage")
    if olgim in [None, "not_mentioned", "0"]:
        return make_result(DecisionStatus.FALSE)

    if not _olgim_low(data):  # 必须 I 或 II
        return make_result(DecisionStatus.FALSE)

    if _olga_high(data):
        return make_result(DecisionStatus.FALSE)

    if has_additional_risk_factor(data):
        return make_result(DecisionStatus.TRUE, extra={"chain": 4})
    return make_result(DecisionStatus.FALSE)


# =============================================================================
# CHAIN 5 – NON‑VISIBLE HIGH‑GRADE DYSPLASIA
# =============================================================================
def evaluate_chain_5(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Non‑visible lesion with HIGH‑GRADE dysplasia.
    Follow‑up in 6 months (REC.21).
    """
    missing = require_fields(data, [
        "dysplasia",
        "dysplasia_grade",
        "dysplasia_visibility",
        "neoplastic_lesion_visibility"
    ])
    if missing:
        return make_result(DecisionStatus.UNKNOWN, missing)

    # 新增：确保 dysplasia 为阳性
    if data.get("dysplasia") not in ["yes", "highly_likely"]:
        return make_result(DecisionStatus.FALSE)

    if data.get("dysplasia_visibility") == "visible":
        return make_result(DecisionStatus.FALSE)
    if data.get("neoplastic_lesion_visibility") == "visible":
        return make_result(DecisionStatus.FALSE)

    if data.get("dysplasia_grade") == "high_grade":
        return make_result(DecisionStatus.TRUE, extra={"chain": 5})
    return make_result(DecisionStatus.FALSE)


# =============================================================================
# CHAIN 6 – NON‑VISIBLE LOW‑GRADE/INDEFINTE DYSPLASIA
# =============================================================================
def evaluate_chain_6(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Non‑visible lesion with LOW‑grade dysplasia OR indefinite for dysplasia
    OR highly_likely dysplasia (without high-grade specified).
    Follow‑up in 12 months (REC.21).
    """
    missing = require_fields(data, [
        "dysplasia",
        "dysplasia_grade",
        "dysplasia_visibility",
        "neoplastic_lesion_visibility"
    ])
    if missing:
        return make_result(DecisionStatus.UNKNOWN, missing)

    if data.get("dysplasia_visibility") == "visible":
        return make_result(DecisionStatus.FALSE)
    if data.get("neoplastic_lesion_visibility") == "visible":
        return make_result(DecisionStatus.FALSE)

    # 1. 若 indefinite_for_dysplasia 为 yes，直接触发链6
    if _yes_no_to_bool(data.get("indefinite_for_dysplasia")) is True:
        return make_result(DecisionStatus.TRUE, extra={"chain": 6})

    # 2. 否则，必须 dysplasia 为阳性
    if data.get("dysplasia") not in ["yes", "highly_likely"]:
        return make_result(DecisionStatus.FALSE)

    cond = (
        data.get("dysplasia_grade") == "low_grade" or
        data.get("dysplasia") == "highly_likely"
    )
    if cond:
        return make_result(DecisionStatus.TRUE, extra={"chain": 6})
    return make_result(DecisionStatus.FALSE)


# =============================================================================
# CHAIN 7 – INDEFINITE / SUSPICIOUS WITH VISIBLE LESION
# =============================================================================
def evaluate_chain_7(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Visible lesion with indefinite for dysplasia OR highly_likely neoplastic/dysplasia.
    Requires expert GI pathology review (REC.22).
    """
    missing = require_fields(data, [
        "indefinite_for_dysplasia",
        "neoplastic_lesion",
        "neoplastic_lesion_visibility",
        "dysplasia",
        "dysplasia_visibility"
    ])
    if missing:
        return make_result(DecisionStatus.UNKNOWN, missing)

    # 条件1：neoplastic_lesion 高度可疑且可见
    cond_neoplastic = (
        data.get("neoplastic_lesion_visibility") == "visible" and
        data.get("neoplastic_lesion") == "highly_likely"
    )

    # 条件2：dysplasia 高度可疑且可见，或 indefinite_for_dysplasia 为 yes 且可见（dysplasia_visibility 可见）
    cond_dysplasia = (
        data.get("dysplasia_visibility") == "visible" and
        (
            data.get("dysplasia") == "highly_likely" or
            _yes_no_to_bool(data.get("indefinite_for_dysplasia")) is True
        )
    )

    if cond_neoplastic or cond_dysplasia:
        return make_result(DecisionStatus.TRUE, extra={"chain": 7})
    return make_result(DecisionStatus.FALSE)


# =============================================================================
# CHAIN 8 – HIGH‑RISK LESION → SURGERY
# =============================================================================
def evaluate_chain_8(data: Dict[str, Any]) -> Dict[str, Any]:
    # 若已切除，不应进入术前手术链
    if data.get("en_bloc_resection") not in [None, "not_mentioned"]:
        return make_result(DecisionStatus.FALSE)

    missing = require_fields(data, [
        "neoplastic_lesion_size", "differentiation_status", "neoplastic_lesion_ulceration"
    ])
    if missing:
        return make_result(DecisionStatus.UNKNOWN, missing)

    size = _parse_measurement(data["neoplastic_lesion_size"])
    depth = data.get("submucosal_invasion_depth")
    diff = data["differentiation_status"]
    ulcer = _yes_no_to_bool(data["neoplastic_lesion_ulceration"]) is True
    sm_inv = data["submucosal_invasion"] == "yes"

    # ---- 1. 明确深部浸润证据 ----
    deep_sm_suspicion = False

    # A. 病理明确深部浸润 (>500µm)
    if sm_inv and depth == "gt_500um":
        deep_sm_suspicion = True

    # B. 当深度未知时，不武断推断；但若有临床征象（溃疡），应提示补充信息
    if sm_inv and ulcer and depth == "not_mentioned":
        # 返回 UNKNOWN，要求 EUS 或进一步检查明确深度
        return make_result(DecisionStatus.UNKNOWN, ["submucosal_invasion_depth"])

    # ---- 2. 明确的手术指征 ----
    cond_undiff = (diff == "undifferentiated" and (ulcer or (size is not None and size > 20)))
    cond_diff_ulcer = (diff == "differentiated" and ulcer and size is not None and size > 30)
    cond_diff_sm_size = (diff == "differentiated" and sm_inv and size is not None and size > 30)

    cond = deep_sm_suspicion or cond_undiff or cond_diff_ulcer or cond_diff_sm_size

    if cond:
        return make_result(DecisionStatus.TRUE, extra={"chain": 8})

    # 其他深度未知但有 SM 浸润的情况（无溃疡、无其他指征）
    if sm_inv and depth == "not_mentioned" and size is None:
        return make_result(DecisionStatus.UNKNOWN, ["neoplastic_lesion_size", "submucosal_invasion_depth"])
    if sm_inv and depth == "not_mentioned":
        return make_result(DecisionStatus.UNKNOWN, ["submucosal_invasion_depth"])

    return make_result(DecisionStatus.FALSE)


# =============================================================================
# CHAIN 9 – CONSERVATIVE MANAGEMENT
# =============================================================================
def evaluate_chain_9(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Pre-therapy conservative management (Figure 8).
    Either explicitly recommended, or clinically inferred from limited life expectancy.
    """
    # 1. Explicit conservative recommendation from report
    rec = data.get("recommended_treatment")
    if rec == "conservative_management":
        return make_result(DecisionStatus.TRUE, extra={"chain": 9, "basis": "explicit_recommendation"})

    # 2. Infer conservative management if patient has limited life expectancy (<10 years)
    #    AND has a visible lesion requiring action, with no prior resection.
    life = data.get("life_expectancy")
    if life == "lt_10_years":
        lesion = data.get("neoplastic_lesion")
        visibility = data.get("neoplastic_lesion_visibility")
        if lesion == "yes" and visibility == "visible":
            if data.get("en_bloc_resection") in [None, "not_mentioned"]:
                return make_result(DecisionStatus.TRUE, extra={"chain": 9, "basis": "limited_life_expectancy"})

    return make_result(DecisionStatus.FALSE)


# =============================================================================
# CHAIN 10 – CURATIVE / VERY LOW-RISK RESECTION
# =============================================================================
def evaluate_chain_10(data: Dict[str, Any]) -> Dict[str, Any]:
    missing = require_fields(data, [
        "en_bloc_resection", "en_bloc_resection_margins", "lymphovascular_invasion",
        "differentiation_status", "submucosal_invasion", "neoplastic_lesion_ulceration",
        "neoplastic_lesion_size"
    ])
    if missing:
        return make_result(DecisionStatus.UNKNOWN, missing)

    size = _parse_measurement(data["neoplastic_lesion_size"])
    ulcer = _yes_no_to_bool(data["neoplastic_lesion_ulceration"]) is True
    sm_inv = data["submucosal_invasion"] == "yes"
    if sm_inv:
        return make_result(DecisionStatus.FALSE)

    base = (_yes_no_to_bool(data["en_bloc_resection"]) is True and
            data["en_bloc_resection_margins"] == "r0" and
            _yes_no_to_bool(data["lymphovascular_invasion"]) is False and
            data["differentiation_status"] == "differentiated")
    if not base:
        return make_result(DecisionStatus.FALSE)

    if ulcer and (size is None or size > 30):
        return make_result(DecisionStatus.FALSE)

    extra = {"chain": 10}
    if emr_eligible(data):
        extra["emr_eligible"] = True
    return make_result(DecisionStatus.TRUE, extra=extra)


# =============================================================================
# CHAIN 11 – CURATIVE / LOW-RISK RESECTION
# =============================================================================
def evaluate_chain_11(data: Dict[str, Any]) -> Dict[str, Any]:
    missing = require_fields(data, [
        "en_bloc_resection", "en_bloc_resection_margins", "lymphovascular_invasion",
        "differentiation_status", "submucosal_invasion", "submucosal_invasion_depth",
        "neoplastic_lesion_size", "neoplastic_lesion_ulceration"
    ])
    if missing:
        return make_result(DecisionStatus.UNKNOWN, missing)

    base = (_yes_no_to_bool(data["en_bloc_resection"]) is True and
            data["en_bloc_resection_margins"] == "r0" and
            _yes_no_to_bool(data["lymphovascular_invasion"]) is False)
    if not base:
        return make_result(DecisionStatus.FALSE)

    size = _parse_measurement(data["neoplastic_lesion_size"])
    depth = data.get("submucosal_invasion_depth")
    diff = data["differentiation_status"]
    ulcer = _yes_no_to_bool(data["neoplastic_lesion_ulceration"]) is True
    sm_inv = data["submucosal_invasion"] == "yes"

    cond_a = False
    if sm_inv and diff == "differentiated" and size is not None and size <= 30:
        if depth == "le_500um":
            cond_a = True
        elif depth == "not_mentioned":
            return make_result(DecisionStatus.UNKNOWN, ["submucosal_invasion_depth"])

    cond_b = (not sm_inv and diff == "undifferentiated" and
              size is not None and size <= 20 and not ulcer)

    if cond_a or cond_b:
        extra = {"chain": 11}
        if emr_eligible(data):
            extra["emr_eligible"] = True
        return make_result(DecisionStatus.TRUE, extra=extra)
    return make_result(DecisionStatus.FALSE)


# =============================================================================
# CHAIN 12 – LOCAL-RISK RESECTION
# =============================================================================
def evaluate_chain_12(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Local-risk resection (Figure 9).
    Piecemeal resection or positive horizontal margin of a lesion otherwise
    meeting very low-risk or low-risk criteria.
    Includes:
    - Very low-risk: differentiated, pT1a, LVI negative, no ulcer or ulcer≤30mm
    - Low-risk (differentiated): pT1b, SM≤500µm, size≤30mm, LVI negative
    - Low-risk (undifferentiated): pT1a, size≤20mm, no ulcer, LVI negative
    """
    missing = require_fields(data, [
        "piecemeal_resection",
        "horizontal_margin_status",
        "differentiation_status",
        "lymphovascular_invasion",
        "submucosal_invasion",
        "neoplastic_lesion_ulceration",
        "neoplastic_lesion_size"
    ])
    if missing:
        return make_result(DecisionStatus.UNKNOWN, missing)

    piecemeal = _yes_no_to_bool(data["piecemeal_resection"]) is True
    hm_positive = data["horizontal_margin_status"] == "positive"
    if not (piecemeal or hm_positive):
        return make_result(DecisionStatus.FALSE)

    # ---- 核心条件：必须无脉管侵犯（LVI阴性） ----
    if _yes_no_to_bool(data["lymphovascular_invasion"]) is not False:
        return make_result(DecisionStatus.FALSE)

    diff = data["differentiation_status"]
    sm_inv = data["submucosal_invasion"] == "yes"
    ulcer = _yes_no_to_bool(data["neoplastic_lesion_ulceration"]) is True
    size = _parse_measurement(data["neoplastic_lesion_size"])
    if size is None:
        return make_result(DecisionStatus.UNKNOWN, ["neoplastic_lesion_size"])

    # ---- 1. 极低风险：分化型，pT1a ----
    if diff == "differentiated" and not sm_inv:
        if (not ulcer) or (ulcer and size <= 30):
            return make_result(DecisionStatus.TRUE, extra={"chain": 12})

    # ---- 2. 低风险（分化型）：pT1b，SM≤500µm，≤30mm ----
    if diff == "differentiated" and sm_inv:
        # 深度必须明确且为≤500µm
        depth = data.get("submucosal_invasion_depth")
        if depth == "le_500um" and size <= 30:
            return make_result(DecisionStatus.TRUE, extra={"chain": 12})
        # 若深度缺失或为gt_500um，则不满足低风险条件
        else:
            return make_result(DecisionStatus.FALSE)

    # ---- 3. 低风险（未分化型）：pT1a，≤20mm，无溃疡 ----
    if diff == "undifferentiated" and not sm_inv and not ulcer and size <= 20:
        return make_result(DecisionStatus.TRUE, extra={"chain": 12})

    # ---- 其他情况均不符合 ----
    return make_result(DecisionStatus.FALSE)


# =============================================================================
# CHAIN 13 – HIGH-RISK (NONCURATIVE) RESECTION
# =============================================================================
def evaluate_chain_13(data: Dict[str, Any]) -> Dict[str, Any]:
    missing = require_fields(data, [
        "vertical_margin_status", "lymphovascular_invasion", "submucosal_invasion_depth",
        "differentiation_status", "neoplastic_lesion_ulceration", "neoplastic_lesion_size",
        "submucosal_invasion"
    ])
    if missing:
        return make_result(DecisionStatus.UNKNOWN, missing)

    depth = data.get("submucosal_invasion_depth")
    size = _parse_measurement(data["neoplastic_lesion_size"])
    diff = data["differentiation_status"]
    ulcer = _yes_no_to_bool(data["neoplastic_lesion_ulceration"]) is True
    sm_inv = data["submucosal_invasion"] == "yes"

    cond_a = (data["vertical_margin_status"] == "positive" or
              _yes_no_to_bool(data["lymphovascular_invasion"]) is True or
              (sm_inv and depth == "gt_500um"))
    # 若深度不明但其他条件已满足，则仍可判定
    if sm_inv and depth == "not_mentioned" and not (data["vertical_margin_status"] == "positive" or
                                                    _yes_no_to_bool(data["lymphovascular_invasion"]) is True):
        return make_result(DecisionStatus.UNKNOWN, ["submucosal_invasion_depth"])

    cond_b = (diff == "undifferentiated" and (ulcer or (size is not None and size > 20)))
    cond_c = False
    if sm_inv and diff == "differentiated" and depth == "le_500um" and size is not None and size > 30:
        cond_c = True
    elif sm_inv and diff == "differentiated" and depth == "not_mentioned" and size is not None and size > 30:
        # 深度不明，无法判断是否为高风险（需区分 ≤500µm vs >500µm）
        # 返回 UNKNOWN 要求补充深度信息，而非直接判定为高风险
        return make_result(DecisionStatus.UNKNOWN, ["submucosal_invasion_depth"])

    cond_d = (not sm_inv and ulcer and size is not None and size > 30)

    if cond_a or cond_b or cond_c or cond_d:
        extra = {"chain": 13}
        if emr_eligible(data):
            extra["emr_eligible"] = True
        return make_result(DecisionStatus.TRUE, extra=extra)
    return make_result(DecisionStatus.FALSE)


# =============================================================================
# CHAIN 14 – AUTOIMMUNE GASTRITIS
# =============================================================================
def evaluate_chain_14(data: Dict[str, Any]) -> Dict[str, Any]:
    missing = require_fields(data, ["autoimmune_gastritis"])
    if missing:
        return make_result(DecisionStatus.UNKNOWN, missing)
    if _yes_no_to_bool(data["autoimmune_gastritis"]) is True:
        return make_result(DecisionStatus.TRUE, extra={"chain": 14})
    return make_result(DecisionStatus.FALSE)


# =============================================================================
# CHAIN 15 – COMMON VARIABLE IMMUNODEFICIENCY
# =============================================================================
def evaluate_chain_15(data: Dict[str, Any]) -> Dict[str, Any]:
    missing = require_fields(data, ["common_variable_immunodeficiency"])
    if missing:
        return make_result(DecisionStatus.UNKNOWN, missing)
    if _yes_no_to_bool(data["common_variable_immunodeficiency"]) is True:
        return make_result(DecisionStatus.TRUE, extra={"chain": 15})
    return make_result(DecisionStatus.FALSE)


# =============================================================================
# CHAIN 16 & 17 – MALT IN REMISSION (动态判断癌前病变)
# =============================================================================
def _has_precancerous_conditions(data: Dict[str, Any]) -> bool:
    """
    判断是否存在需要监测的癌前病变（即高危/广泛改变）。
    轻度萎缩（K-T C0-2、OLGA 0-II）视为无癌前病变。
    """
    # 1. 如果存在高危特征，则为有癌前病变
    if _kt_high(data) or _eggim_high(data) or _olga_high(data) or _olgim_high(data):
        return True
    if _yes_no_to_bool(data.get("incomplete_intestinal_metaplasia")) is True:
        return True

    # 2. 如果所有客观指标均为明确的低风险（或未提及），视为无癌前病变
    kt_low = _kt_low(data)
    olga_low = _olga_low(data)
    olgim_low = _olgim_low(data)
    eggim_low = _eggim_low(data)
    incomplete = _yes_no_to_bool(data.get("incomplete_intestinal_metaplasia"))

    # 只要低风险指标（含未提及）且无高危，返回 False
    if (kt_low is not False and olga_low is not False and
        olgim_low is not False and eggim_low is not False and
        incomplete is not True):
        return False

    # 其他情况：有明确高危或信息不足，保守起见视为有癌前病变
    return True

def evaluate_chain_16(data: Dict[str, Any]) -> Dict[str, Any]:
    missing = require_fields(data, ["gastric_malt_lymphoma"])
    if missing:
        return make_result(DecisionStatus.UNKNOWN, missing)
    if data["gastric_malt_lymphoma"] == "in_remission" and not _has_precancerous_conditions(data):
        return make_result(DecisionStatus.TRUE, extra={"chain": 16})
    return make_result(DecisionStatus.FALSE)

def evaluate_chain_17(data: Dict[str, Any]) -> Dict[str, Any]:
    missing = require_fields(data, ["gastric_malt_lymphoma"])
    if missing:
        return make_result(DecisionStatus.UNKNOWN, missing)
    if data["gastric_malt_lymphoma"] == "in_remission" and _has_precancerous_conditions(data):
        return make_result(DecisionStatus.TRUE, extra={"chain": 17})
    return make_result(DecisionStatus.FALSE)


# =============================================================================
# CHAIN 18 – HEREDITARY SYNDROME
# =============================================================================
def evaluate_chain_18(data: Dict[str, Any]) -> Dict[str, Any]:
    missing = require_fields(data, ["hereditary_syndrome"])
    if missing:
        return make_result(DecisionStatus.UNKNOWN, missing)
    if _yes_no_to_bool(data["hereditary_syndrome"]) is True:
        return make_result(DecisionStatus.TRUE, extra={"chain": 18})
    return make_result(DecisionStatus.FALSE)


# =============================================================================
# CHAIN 19 – SEROLOGIC HIGH RISK
# =============================================================================
def evaluate_chain_19(data: Dict[str, Any]) -> Dict[str, Any]:
    missing = require_fields(data, ["pepsinogen_i_level", "pepsinogen_ii_level"])
    if missing:
        return make_result(DecisionStatus.UNKNOWN, missing)
    if high_serologic_risk(data):
        return make_result(DecisionStatus.TRUE, extra={"chain": 19})
    return make_result(DecisionStatus.FALSE)


# =============================================================================
# CHAIN 20 – STOP SURVEILLANCE (AGE / LIFE EXPECTANCY)
# =============================================================================
def evaluate_chain_20(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Guideline 6 & 23: Stop surveillance in ASYMPTOMATIC individuals over 80
    or with life expectancy <10 years.
    If there is a visible lesion or dysplasia, do NOT stop.
    """
    # 若存在可见的肿瘤或异型增生，则不应停止监测
    lesion = data.get("neoplastic_lesion")
    dys = data.get("dysplasia")
    if lesion in ["yes", "highly_likely"] or dys in ["yes", "highly_likely"]:
        return make_result(DecisionStatus.FALSE)

    age = data.get("patient_age")
    life = data.get("life_expectancy")

    # 如果两个字段都缺失，则无法判断，返回 FALSE
    if (age in [None, "not_mentioned"] and life in [None, "not_mentioned"]):
        return make_result(DecisionStatus.FALSE)

    age_ok = False
    if age not in [None, "not_mentioned"]:
        age_val = _safe_int(age)
        if age_val is not None and age_val >= 80:
            age_ok = True

    life_ok = (life == "lt_10_years")

    if age_ok or life_ok:
        return make_result(DecisionStatus.TRUE, extra={"chain": 20})
    return make_result(DecisionStatus.FALSE)
    
# =============================================================================
# New: Chain 21 (visible lesions suitable for ESD)
# =============================================================================
def evaluate_chain_21(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Pre-therapy visible lesion suitability for ESD (Figure 8).
    Covers both:
    - Dysplastic lesion, any size (if no deep SM invasion)
    - Differentiated/Undifferentiated carcinoma (with size/ulcer limits)
    """
    # 1. Basic conditions: visible, confirmed neoplastic OR dysplastic, not yet resected
    if data.get("en_bloc_resection") not in [None, "not_mentioned"]:
        return make_result(DecisionStatus.FALSE)

    # ----- 检查是否存在可见的异型增生或肿瘤 -----
    is_visible = False
    lesion_type = None  # 用于记录触发的类型

    # 检查 neoplastic_lesion（肿瘤）
    if data.get("neoplastic_lesion") == "yes" and data.get("neoplastic_lesion_visibility") == "visible":
        is_visible = True
        lesion_type = "neoplastic"
    # 检查 dysplasia（异型增生）—— 这是新增的核心部分
    if data.get("dysplasia") == "yes" and data.get("dysplasia_visibility") == "visible":
        is_visible = True
        lesion_type = "dysplastic"

    if not is_visible:
        return make_result(DecisionStatus.FALSE)

    # 2. Mandatory fields: size and ulceration (clinical decision basis)
    missing = require_fields(data, [
        "neoplastic_lesion_size",
        "neoplastic_lesion_ulceration"
    ])
    if missing:
        return make_result(DecisionStatus.UNKNOWN, missing)

    size = _parse_measurement(data["neoplastic_lesion_size"])
    ulcer = _yes_no_to_bool(data["neoplastic_lesion_ulceration"]) is True
    if size is None:
        return make_result(DecisionStatus.UNKNOWN, ["neoplastic_lesion_size"])

    # 3. Get optional fields (may be "not_mentioned")
    diff = data.get("differentiation_status")
    sm_inv = data.get("submucosal_invasion")
    depth = data.get("submucosal_invasion_depth")

    # 4. Exclude if there is clear evidence of deep submucosal invasion
    if sm_inv == "yes" and depth == "gt_500um":
        return make_result(DecisionStatus.FALSE)

    # 5. 如果病变类型为异型增生（dysplasia），直接按分化型黏膜内癌处理（无溃疡大小不限，有溃疡≤30mm）
    #    因为异型增生通常为分化型且未浸润。
    if lesion_type == "dysplastic":
        # 异型增生：无溃疡 -> 大小不限；有溃疡 -> ≤30mm
        if (not ulcer) or (ulcer and size <= 30):
            # 注意：异型增生通常不适用 EMR（除非满足巴黎分型），但仍可标记 emr_eligible
            extra = {"chain": 21, "basis": "dysplastic"}
            if emr_eligible(data):
                extra["emr_eligible"] = True
            return make_result(DecisionStatus.TRUE, extra=extra)
        else:
            return make_result(DecisionStatus.FALSE)

    # 6. 如果病变类型为肿瘤（neoplastic），按原有分化逻辑判断
    # ------------------------------------------------------------
    if diff == "differentiated":
        # 分化型：无溃疡 -> 任何大小；有溃疡 -> ≤30mm
        if (not ulcer) or (ulcer and size <= 30):
            extra = {"chain": 21, "basis": "differentiated"}
            if emr_eligible(data):
                extra["emr_eligible"] = True
            return make_result(DecisionStatus.TRUE, extra=extra)
        else:
            return make_result(DecisionStatus.FALSE)

    elif diff == "undifferentiated":
        # 未分化型：必须黏膜内（sm_inv == "no"），无溃疡，≤20mm
        if sm_inv == "no" and (not ulcer) and size <= 20:
            return make_result(DecisionStatus.TRUE, extra={"chain": 21, "basis": "undifferentiated"})
        else:
            return make_result(DecisionStatus.FALSE)

    else:  # diff == "not_mentioned" or None
        # 分化未知：仅安全情况是无溃疡且≤20mm
        if (not ulcer) and size <= 20:
            extra = {"chain": 21, "basis": "differentiation_unknown_safe"}
            if emr_eligible(data):
                extra["emr_eligible"] = True
            return make_result(DecisionStatus.TRUE, extra=extra)
        # 否则暂不推荐ESD，转手术/MDT评估
        return make_result(DecisionStatus.FALSE)


# =============================================================================
# PRECEDENCE HIERARCHY
# =============================================================================
CHAIN_GROUPS = [
    [20],
    [10,11,12,13,21,8,9],
    [5,6,7],
    [14,16],
    [15,17,18,19],
    [3,2,4,1]
]

CHAIN_FUNCTIONS = {
    1: evaluate_chain_1, 2: evaluate_chain_2, 3: evaluate_chain_3, 4: evaluate_chain_4,
    5: evaluate_chain_5, 6: evaluate_chain_6, 7: evaluate_chain_7, 8: evaluate_chain_8,
    9: evaluate_chain_9, 10: evaluate_chain_10, 11: evaluate_chain_11, 12: evaluate_chain_12,
    13: evaluate_chain_13, 14: evaluate_chain_14, 15: evaluate_chain_15, 16: evaluate_chain_16,
    17: evaluate_chain_17, 18: evaluate_chain_18, 19: evaluate_chain_19, 20: evaluate_chain_20, 21: evaluate_chain_21
}


# =============================================================================
# 并行校验函数（pt_category, resection_curability, recommended_treatment, LNM, atrophic）
# =============================================================================
def parse_lnm_risk(data: Dict[str, Any]) -> Optional[str]:
    risk = data.get("lymph_node_metastasis_risk")
    if risk is None or risk == "not_mentioned":
        return None
    if risk in ["lt_0.5_percent", "0.5_to_1_percent"]:
        return "very_low"
    elif risk == "lt_3_percent":
        return "low"
    elif risk == "ge_3_percent":
        return "high"
    return None

def validate_pt_category(data: Dict[str, Any], chain_result: Dict[str, Any]) -> List[str]:
    warnings = []
    pt = data.get("pt_category")
    if pt is None or pt == "not_mentioned":
        return warnings
    sm = data.get("submucosal_invasion")
    depth = data.get("submucosal_invasion_depth")
    if pt == "pt1a" and sm == "yes":
        warnings.append("pt_category is pt1a but submucosal_invasion is yes, inconsistent.")
    elif pt == "pt1b" and sm == "no":
        warnings.append("pt_category is pt1b but submucosal_invasion is no, inconsistent.")
    if pt == "pt1b" and depth == "not_mentioned":
        warnings.append("pt_category is pt1b but submucosal_invasion_depth is not mentioned, please clarify.")
    return warnings

def validate_resection_curability(data: Dict[str, Any], chain_result: Dict[str, Any]) -> List[str]:
    warnings = []
    curability = data.get("resection_curability")
    if curability is None or curability == "not_mentioned":
        return warnings
    chain = chain_result.get("chain")
    if chain is None:
        return warnings
    if chain in [10, 11] and curability == "noncurative":
        warnings.append(f"Chain {chain} indicates curative resection, but resection_curability is noncurative, mismatch.")
    elif chain == 13 and curability == "curative":
        warnings.append(f"Chain {chain} indicates noncurative resection, but resection_curability is curative, mismatch.")
    return warnings

def validate_recommended_treatment(data: Dict[str, Any], chain_result: Dict[str, Any]) -> List[str]:
    """
    校验 recommended_treatment 与指南标准治疗是否一致。
    仅对治疗前链（8,9,21）进行校验，术后链（10-13）不校验。
    允许 'endoscopic_resection' 作为 'endoscopic_submucosal_dissection' 或 'endoscopic_mucosal_resection' 的等价表示。
    """
    warnings = []
    rec = data.get("recommended_treatment")
    if rec is None or rec == "not_mentioned":
        return warnings
    
    chain = chain_result.get("chain")
    if chain is None:
        return warnings
    
    # 仅对治疗前链（8,9,21）进行校验
    if chain not in [8, 9, 21]:
        return warnings    
    # 确定该链对应的指南标准治疗
    if chain == 8:
        standard = "surgical_treatment"
    elif chain == 21:
        # 如果 emr_eligible 为 True，标准治疗可以是 EMR（替代ESD）
        if chain_result.get("extra", {}).get("emr_eligible"):
            standard = "endoscopic_mucosal_resection"
        else:
            standard = "endoscopic_submucosal_dissection"
    else:  # chain == 9
        standard = "conservative_management"    
   
    # 检查是否一致：允许 endoscopic_resection 作为 ESD/EMR 的同义词
    if chain == 21:
        # 如果 rec 是 "endoscopic_resection"，视为与任何内镜切除标准一致
        if rec == "endoscopic_resection":
            return warnings  # 不报错
        # 否则比较 exact match
        if rec != standard:
            warnings.append(
                f"Suggested treatment '{rec}' differs from guideline standard '{standard}', MDT review advised."
            )
    else:
        # 其他链直接比较
        if rec not in [standard, "not_mentioned", "conservative_management"]:
            warnings.append(
                f"Suggested treatment '{rec}' differs from guideline standard '{standard}', MDT review advised."
            )
    
    return warnings

def validate_lnm_risk(data: Dict[str, Any], chain_result: Dict[str, Any]) -> List[str]:
    warnings = []
    lnm = data.get("lymph_node_metastasis_risk")
    if lnm is None or lnm == "not_mentioned":
        return warnings
    lnm_level = parse_lnm_risk(data)
    chain = chain_result.get("chain")
    if chain is None:
        return warnings
    if chain == 10:
        chain_level = "very_low"
    elif chain == 11:
        chain_level = "low"
    elif chain == 13:
        chain_level = "high"
    elif chain == 12:
        return warnings  # 不比较
    else:
        return warnings
    if lnm_level and lnm_level != chain_level:
        warnings.append(f"Chain {chain} risk level '{chain_level}' conflicts with reported LNM risk '{lnm}', please verify.")
    return warnings

def validate_atrophic_gastritis(data: Dict[str, Any]) -> List[str]:
    warnings = []
    atrophic = data.get("atrophic_gastritis")
    if atrophic is None or atrophic == "not_mentioned":
        return warnings
    kt = data.get("kimura_takemoto_classification")
    olga = data.get("olga_stage")
    high = False
    if kt and kt in ["c3","o1","o2","o3"]:
        high = True
    if olga and olga in ["iii","iv"]:
        high = True
    if high and atrophic == "no":
        warnings.append("atrophic_gastritis is 'no' but Kimura/OLGA suggests significant atrophy, please review.")
    return warnings


# =============================================================================
# 主分类器
# =============================================================================
# MODIFIER CHAINS (这些链不互斥，应叠加在基础链之上)
MODIFIER_CHAINS = {15, 17, 18, 19}
BASE_CHAINS = {1, 2, 3, 4}

# 内部辅助：生成治疗建议（仅针对 Chain 21 的 EMR 替代方案）
def _get_treatment_suggestion(chain, details):
    if chain is None:
        return None
    
    # 处理组合链（如 "21,2" 或 "17,21" 等）
    chain_str = str(chain)
    if isinstance(chain, str) and "," in chain_str:
        chain_ids = [int(c.strip()) for c in chain_str.split(",") if c.strip().isdigit()]
    else:
        chain_ids = [chain] if isinstance(chain, int) else []
    
    # 如果链中包含 21，检查 emr_eligible 标记
    if 21 in chain_ids:
        # 尝试从 details 中提取 extra
        extra = {}
        if isinstance(details, dict):
            # 组合链时 details 结构为 {"modifier": {...}, "base": {...}, "extra": {...}}
            extra = details.get("extra", {})
            if not extra:
                # 非组合链时 details 直接包含 extra
                extra = details.get("extra", {})
        if extra.get("emr_eligible"):
            return "EMR is an alternative for Paris 0-IIa lesions ≤10mm with low malignancy likelihood.(REC.26)"
    return None
    
def classify_patient(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    主分类器，逻辑：
    1. 收集所有链的 TRUE 和 UNKNOWN。
    2. 如果有治疗/术后链（10,11,12,13,21,8,9），优先返回（可附加修饰链）。
    3. 否则，如果有基础链（1-4），智能选择（优先链3）。
    4. 可附加多个修饰链（15,17,18,19）。
    5. 安全网：仅当有可见病变且所有治疗链均无法判定（UNKNOWN）时，才提示缺失字段，而非直接 N/A。
    """
    unknown_fields: List[str] = []
    all_positives: List[Dict[str, Any]] = []
    modifier_positives: List[Dict[str, Any]] = []
    base_positives: List[Dict[str, Any]] = []
    treatment_positives: List[Dict[str, Any]] = []
    dysplasia_positives: List[Dict[str, Any]] = []
    special_positives: List[Dict[str, Any]] = []
    stop_positives: List[Dict[str, Any]] = []

    # 1. 遍历所有链，收集 TRUE 和 UNKNOWN
    for group in CHAIN_GROUPS:
        for chain_id in group:
            result = CHAIN_FUNCTIONS[chain_id](data)
            if result["status"] == DecisionStatus.UNKNOWN:
                unknown_fields.extend(result.get("missing_fields", []))
            elif result["status"] == DecisionStatus.TRUE:
                all_positives.append(result)
                if chain_id in MODIFIER_CHAINS:
                    modifier_positives.append(result)
                elif chain_id in BASE_CHAINS:
                    base_positives.append(result)
                elif chain_id in {21, 8, 9, 10, 11, 12, 13}:
                    treatment_positives.append(result)
                elif chain_id in {5, 6, 7}:
                    dysplasia_positives.append(result)
                elif chain_id in {14, 16}:
                    special_positives.append(result)
                elif chain_id == 20:
                    stop_positives.append(result)

    # 2. 治疗/术后链优先（不与其他基础链组合，但可附加修饰链）
    if treatment_positives:
        # 选取第一个治疗链（按优先级，通常只有一个）
        treatment = treatment_positives[0]
        # 检查是否有修饰链可附加
        modifiers = [p for p in all_positives if p["chain"] in MODIFIER_CHAINS]
        # 构建输出
        chain = treatment["chain"]
        chain_list = [chain] + [m["chain"] for m in modifiers]
        combined_chain = ",".join(str(c) for c in chain_list)

        # 生成建议文本
        tx_text = get_chain_result(chain)
        mod_texts = [get_chain_result(m["chain"]) for m in modifiers]
        surv_text = " | ".join(mod_texts + [tx_text]) if mod_texts else tx_text

        # 合并警告
        warnings = []
        for res in [treatment] + modifiers:
            warnings.extend(validate_pt_category(data, res))
            warnings.extend(validate_resection_curability(data, res))
            warnings.extend(validate_recommended_treatment(data, res))
            warnings.extend(validate_lnm_risk(data, res))
        warnings.extend(validate_atrophic_gastritis(data))

        # 合并 extra
        extra = {}
        for res in [treatment] + modifiers:
            extra.update(res.get("extra", {}))
        extra["combined"] = True

        return {
            "status": "DETERMINISTIC",
            "chain": combined_chain,
            "chain_list": chain_list,
            "details": {"treatment": treatment, "modifiers": modifiers, "extra": extra},
            "missing_fields": sorted(list(set(unknown_fields))),
            "warnings": warnings,
            "surveillance_text": surv_text,
            "treatment_suggestion": _get_treatment_suggestion(combined_chain,
                                                             {"treatment": treatment, "modifiers": modifiers, "extra": extra}),
            "hp_eradication_suggested": (data.get("helicobacter_pylori_status") == "positive"),
            "smoking_cessation_suggested": False,
            "advanced_olg_staging": treatment.get("extra", {}).get("advanced_olg_staging", False)
        }
    # 2.5 停止监测链（链20）优先于异型增生链和基础链
    if stop_positives:
        stop = stop_positives[0]  # 通常只有一个
        return {
            "status": "DETERMINISTIC",
            "chain": 20,
            "chain_list": [20],
            "details": stop,
            "missing_fields": sorted(list(set(unknown_fields))),
            "warnings": [],
            "surveillance_text": get_chain_result(20),
            "treatment_suggestion": None,
            "hp_eradication_suggested": False,
            "smoking_cessation_suggested": False,
            "advanced_olg_staging": False
        }

    # 2.6 术后安全网（在异型增生链之前）
    if data.get("en_bloc_resection") not in [None, "not_mentioned"]:
        postoperative_true = any(p["chain"] in [10,11,12,13] for p in all_positives)
        if not postoperative_true:
            missing_postop = []
            for f in ["submucosal_invasion", "submucosal_invasion_depth", "lymphovascular_invasion", "vertical_margin_status"]:
                if data.get(f) in [None, "not_mentioned"]:
                    missing_postop.append(f)
            return {
                "status": "N/A",
                "missing_fields": sorted(list(set(unknown_fields + missing_postop))),
                "warnings": ["Patient has undergone resection but postoperative pathological details are insufficient for risk stratification (Chains 10-13). Please provide complete resection pathology."],
                "chain": None,
                "chain_list": [],
                "details": None,
                "surveillance_text": None,
                "treatment_suggestion": None,
                "hp_eradication_suggested": False,
                "smoking_cessation_suggested": False,
                "advanced_olg_staging": False
            }


    # 3. 异型增生链（5,6,7）优先于特殊链和基础链
    if dysplasia_positives and not treatment_positives:
        chain_list = [p["chain"] for p in dysplasia_positives]
        combined_chain = ",".join(str(c) for c in chain_list)
        surv_text = " | ".join(get_chain_result(c) for c in chain_list)

        warnings = []
        for res in dysplasia_positives:
            warnings.extend(validate_pt_category(data, res))
            warnings.extend(validate_resection_curability(data, res))
            warnings.extend(validate_recommended_treatment(data, res))
            warnings.extend(validate_lnm_risk(data, res))
        warnings.extend(validate_atrophic_gastritis(data))

        extra = {}
        for res in dysplasia_positives:
            extra.update(res.get("extra", {}))
        extra["combined"] = True

        return {
            "status": "DETERMINISTIC",
            "chain": combined_chain,
            "chain_list": chain_list,
            "details": {"dysplasia": dysplasia_positives, "extra": extra},
            "missing_fields": sorted(list(set(unknown_fields))),
            "warnings": warnings,
            "surveillance_text": surv_text,
            "treatment_suggestion": _get_treatment_suggestion(combined_chain, {"dysplasia": dysplasia_positives, "extra": extra}),
            "hp_eradication_suggested": (data.get("helicobacter_pylori_status") == "positive"),
            "smoking_cessation_suggested": False,
            "advanced_olg_staging": False
        }
    
    # 4. 仅有基础链（无治疗链）
    if base_positives and not treatment_positives:
        # 智能选择：优先链3，其次2，4，1
        base = None
        for chain_id in [3, 2, 4, 1]:
            for p in base_positives:
                if p["chain"] == chain_id:
                    base = p
                    break
            if base:
                break
        if not base:
            base = base_positives[0]

        chain = base["chain"]
        
        # ===== 新增：如果选中的基础链不是链3，但存在特殊链，则转用特殊链 =====
        if chain != 3 and special_positives:
            # 直接复用特殊链处理逻辑（附加修饰链）
            modifiers = [p for p in all_positives if p["chain"] in MODIFIER_CHAINS]
            special_chain_list = [p["chain"] for p in special_positives] + [m["chain"] for m in modifiers]
            combined_chain = ",".join(str(c) for c in special_chain_list)
            surv_text = " | ".join([get_chain_result(c) for c in special_chain_list])
            
            warnings = []
            for res in special_positives + modifiers:
                warnings.extend(validate_pt_category(data, res))
                warnings.extend(validate_resection_curability(data, res))
                warnings.extend(validate_recommended_treatment(data, res))
                warnings.extend(validate_lnm_risk(data, res))
            warnings.extend(validate_atrophic_gastritis(data))
            
            extra = {}
            for res in special_positives + modifiers:
                extra.update(res.get("extra", {}))
            extra["combined"] = True
            
            return {
                "status": "DETERMINISTIC",
                "chain": combined_chain,
                "chain_list": special_chain_list,
                "details": {"special": special_positives, "modifiers": modifiers, "extra": extra},
                "missing_fields": sorted(list(set(unknown_fields))),
                "warnings": warnings,
                "surveillance_text": surv_text,
                "treatment_suggestion": _get_treatment_suggestion(combined_chain, {"special": special_positives, "modifiers": modifiers, "extra": extra}),
                "hp_eradication_suggested": (data.get("helicobacter_pylori_status") == "positive"),
                "smoking_cessation_suggested": False,
                "advanced_olg_staging": False
            }
            
        # 若为链3，或没有特殊链，则继续原有的基础链逻辑（附加修饰链）
        modifiers = [p for p in all_positives if p["chain"] in MODIFIER_CHAINS]
        chain_list = [chain] + [m["chain"] for m in modifiers]
        combined_chain = ",".join(str(c) for c in chain_list)

        tx_text = get_chain_result(chain)
        mod_texts = [get_chain_result(m["chain"]) for m in modifiers]
        surv_text = " | ".join(mod_texts + [tx_text]) if mod_texts else tx_text

        warnings = []
        for res in [base] + modifiers:
            warnings.extend(validate_pt_category(data, res))
            warnings.extend(validate_resection_curability(data, res))
            warnings.extend(validate_recommended_treatment(data, res))
            warnings.extend(validate_lnm_risk(data, res))
        warnings.extend(validate_atrophic_gastritis(data))

        extra = {}
        for res in [base] + modifiers:
            extra.update(res.get("extra", {}))
        extra["combined"] = True

        return {
            "status": "DETERMINISTIC",
            "chain": combined_chain,
            "chain_list": chain_list,
            "details": {"base": base, "modifiers": modifiers, "extra": extra},
            "missing_fields": sorted(list(set(unknown_fields))),
            "warnings": warnings,
            "surveillance_text": surv_text,
            "treatment_suggestion": _get_treatment_suggestion(combined_chain,
                                                               {"base": base, "modifiers": modifiers, "extra": extra}),
            "hp_eradication_suggested": (data.get("helicobacter_pylori_status") == "positive"),
            "smoking_cessation_suggested": False,
            "advanced_olg_staging": base.get("extra", {}).get("advanced_olg_staging", False)
        }

    # 5. 特殊情境链（链14、16）优先于修饰链和基础链，但低于治疗链
    if special_positives and not treatment_positives and not base_positives:
        # 检查是否有修饰链可附加
        modifiers = [p for p in all_positives if p["chain"] in MODIFIER_CHAINS]
        chain_list = [p["chain"] for p in special_positives] + [m["chain"] for m in modifiers]
        combined_chain = ",".join(str(c) for c in chain_list)
        surv_text = " | ".join([get_chain_result(c) for c in chain_list])

        warnings = []
        for res in special_positives + modifiers:
            warnings.extend(validate_pt_category(data, res))
            warnings.extend(validate_resection_curability(data, res))
            warnings.extend(validate_recommended_treatment(data, res))
            warnings.extend(validate_lnm_risk(data, res))
        warnings.extend(validate_atrophic_gastritis(data))

        extra = {}
        for res in special_positives + modifiers:
            extra.update(res.get("extra", {}))
        extra["combined"] = True
    
        return {
            "status": "DETERMINISTIC",
            "chain": combined_chain,
            "chain_list": chain_list,
            "details": {"special": special_positives, "modifiers": modifiers, "extra": extra},
            "missing_fields": sorted(list(set(unknown_fields))),
            "warnings": warnings,
            "surveillance_text": surv_text,
            "treatment_suggestion": _get_treatment_suggestion(combined_chain, {"special": special_positives, "modifiers": modifiers, "extra": extra}),
            "hp_eradication_suggested": (data.get("helicobacter_pylori_status") == "positive"),
            "smoking_cessation_suggested": False,
            "advanced_olg_staging": False
        }

    # 6. 仅有修饰链（无基础链，无治疗链）
    if modifier_positives and not base_positives and not treatment_positives:
        # 多个修饰链可合并
        chain_list = [p["chain"] for p in modifier_positives]
        combined_chain = ",".join(str(c) for c in chain_list)
        surv_text = " | ".join(get_chain_result(c) for c in chain_list)

        warnings = []
        for res in modifier_positives:
            warnings.extend(validate_pt_category(data, res))
            warnings.extend(validate_resection_curability(data, res))
            warnings.extend(validate_recommended_treatment(data, res))
            warnings.extend(validate_lnm_risk(data, res))
        warnings.extend(validate_atrophic_gastritis(data))
        if 17 in chain_list or 18 in chain_list or 19 in chain_list:
            warnings.append(
                "Modifier chain(s) detected but no baseline precancerous staging (Chains 1-4) could be determined. "
                "Please ensure OLGA/OLGIM/Kimura/EGGIM data are available for optimal surveillance interval."
            )

        extra = {}
        for res in modifier_positives:
            extra.update(res.get("extra", {}))
        extra["combined"] = True

        return {
            "status": "DETERMINISTIC",
            "chain": combined_chain,
            "chain_list": chain_list,
            "details": {"modifiers": modifier_positives, "extra": extra},
            "missing_fields": sorted(list(set(unknown_fields))),
            "warnings": warnings,
            "surveillance_text": surv_text,
            "treatment_suggestion": _get_treatment_suggestion(combined_chain, {"modifiers": modifier_positives, "extra": extra}),
            "hp_eradication_suggested": (data.get("helicobacter_pylori_status") == "positive"),
            "smoking_cessation_suggested": False,
            "advanced_olg_staging": any(m.get("extra", {}).get("advanced_olg_staging", False) for m in modifier_positives)
        }

    # 7. AMBIGUOUS 或 N/A
    if len(all_positives) > 1:
        return {
            "status": "AMBIGUOUS",
            "eligible_chains": [x["chain"] for x in all_positives],
            "warnings": [],
            "missing_fields": sorted(list(set(unknown_fields))),
            "chain": None,
            "chain_list": [],
            "details": None,
            "surveillance_text": None,
            "treatment_suggestion": None,
            "hp_eradication_suggested": False,
            "smoking_cessation_suggested": False,
            "advanced_olg_staging": False
        }

    # 真正的 N/A
    return {
        "status": "N/A",
        "missing_fields": sorted(list(set(unknown_fields))),
        "warnings": [],
        "chain": None,
        "chain_list": [],
        "details": None,
        "surveillance_text": None,
        "treatment_suggestion": None,
        "hp_eradication_suggested": False,
        "smoking_cessation_suggested": False,
        "advanced_olg_staging": False
    }


# =============================================================================
# 链结果文本映射
# =============================================================================
def get_chain_result(chain: int) -> str:
    mapping = {
        1: "No endoscopic surveillance recommended after successful H. pylori eradication (if initially positive). (REC.34)",
        2: "Surveillance every 3 years with high-quality endoscopy.(REC.31)",
        3: "Intensive surveillance every 1-2 years with high-quality endoscopy.(REC.33)",
        4: "High quality endoscopic surveillance every 3 years may be considered.(REC.35)",
        5: "Repeat high-quality endoscopy with VCE, staging of precancerous conditions, and H. pylori testing. If again no neoplastic lesion is seen: Follow-up high-quality endoscopy in 6 months (high-grade dysplasia).(REC.21)",
        6: "Perform a high-quality endoscopic re-evaluation. If no lesions are detected: Follow-up endoscopy in 12 months (low-grade dysplasia or indefinite for dysplasia).(REC.21)",
        7: "Referral for expert GI pathologist; high-quality endoscopy with VCE; decision for targeted biopsies or resection based on findings.(REC.22)",
        8: "Surgical treatment (gastrectomy + lymphadenectomy) is recommended. (Fig.8,MAPSIII,REC.24)",
        9: "MDT discussion, consider conservative management.(Fig.8,MAPSIII,REC.23)",
        10: "No further staging procedure or treatment is recommended (REC.28). Surveillance endoscopy at 3-6 months, then annually; if metachronous (or synchronous) lesion detected, similar approach as for any primary gastric lesion (REC.29, Fig.9). Routine imaging (EUS/CT/MRI/PET) is NOT suggested for follow-up of very low-risk resections. (REC.29)",
        11: "Staging should be completed, and further treatment is generally not necessary after a multidisciplinary discussion (REC.28). Surveillance endoscopy at 3-6 months, then annually; if metachronous (or synchronous) lesion detected, similar approach as for any primary gastric lesion (REC.29, Fig.9). ",
        12: "Endoscopic surveillance/re-treatment (repeat ESD, argon plasma coagulation [APC], surgery) is recommended rather than other additional treatment, surveillance endoscopy at 3-6 months, then annually if no recurrence(REC.28-29, Fig.9).",
        13: "Complete staging and multidisciplinary discussion. Consider surgery according to patient characteristics (taking into account age, comorbidities, and life expectancy) and consider LNM risk according to eCura2/W-eCura score(REC.28, Fig.9).if Unfit for patient, consider conservative management.(REC.30)",
        14: "Autoimmune gastritis. High quality endoscopic follow-up every 3 years to detect gastric cancer and neuroendocrine tumors.(REC.47)",
        15: "CVID. Perform high-quality endoscopy at diagnosis; subsequent follow-up based on stage of precancerous conditions and/or autoimmune gastritis.(REC.48)",
        16: "MALT lymphoma in remission, no precancerous conditions. Surveillance every 5 years.(text: expert opinion, p540)",
        17: "MALT lymphoma in remission, with precancerous conditions. Followed up according to the stage of precancerous conditions.(text: expert opinion, p540)",
        18: "Hereditary syndrome. Endoscopic surveillance should follow recommendations for the specific syndrome or according to the gastric mucosal changes, whichever interval is shorter.(REC.46)",
        19: "Individuals with Low Serum Pepsinogen. Endoscopic screening for precancerous conditions recommended.(REC.7)",
        20: "Discontinue or do not start gastric cancer screening or surveillance.(REC.6, REC.23)",
        21: "ESD is recommended based on dysplasia/carcinoma.(REC.24-27)"
    }
    return mapping.get(chain, "Unknown")


# =============================================================================
# 最终输出格式化函数（供 integrate 调用）
# =============================================================================
def format_final_output(classification_result: Dict[str, Any], data: Dict[str, Any]) -> str:
    lines = []
    lines.append("MAPS III 2025 Gastric Cancer Surveillance Decision Report:")

    status = classification_result.get("status")
    if status == "AMBIGUOUS":
        lines.append("Status: AMBIGUOUS – multiple chains satisfied, MDT review required.")
        lines.append(f"Eligible chains: {', '.join(map(str, classification_result.get('eligible_chains', [])))}")
        if classification_result.get("missing_fields"):
            lines.append(f"Missing fields: {', '.join(classification_result['missing_fields'])}")
        return "\n".join(lines)

    if status == "N/A":
        lines.append("Status: N/A – insufficient key information.")
        if classification_result.get("missing_fields"):
            lines.append(f"Missing fields: {', '.join(classification_result['missing_fields'])}")
        return "\n".join(lines)

    chain = classification_result.get("chain")
    lines.append(f"Status: DETERMINISTIC")
    lines.append(f"Chain: {chain}")
    # 优先使用组合建议文本
    surv_text = classification_result.get("surveillance_text")
    if surv_text:
        lines.append(f"Surveillance recommendation: {surv_text}")
    elif chain is not None:
        # 若 chain 是组合字符串（如 "17,2"），需分别获取并拼接
        if isinstance(chain, str) and "," in chain:
            parts = chain.split(",")
            texts = [get_chain_result(int(p.strip())) for p in parts if p.strip().isdigit()]
            lines.append(f"Surveillance recommendation: {' | '.join(texts)}")
        else:
            lines.append(f"Surveillance recommendation: {get_chain_result(chain)}")
    # 动态 Guideline 37 提示（仅当基础链为2或3时）
    # 注意：组合链时 advanced_olg_staging 保存在 base 中
    chain_list = classification_result.get("chain_list", [])
    if classification_result.get("advanced_olg_staging") and (2 in chain_list or 3 in chain_list):
        lines.append("\n Note: This patient has advanced OLGA/OLGIM stages (III/IV). During surveillance, random biopsies are not required if no visible lesions are observed (REC.37).")

    # Decision basis (key field snapshot)
    key_fields = [
        "kimura_takemoto_classification", "eggim_score", "olga_stage", "olgim_stage",
        "incomplete_intestinal_metaplasia", "family_history", "helicobacter_pylori_status",
        "dysplasia", "dysplasia_grade", "neoplastic_lesion_size", "differentiation_status",
        "submucosal_invasion", "submucosal_invasion_depth", "lymphovascular_invasion",
        "en_bloc_resection", "horizontal_margin_status", "vertical_margin_status"
    ]
    lines.append("\n Decision basis (key field values):")
    for f in key_fields:
        val = data.get(f)
        if val is not None and val != "not_mentioned":
            lines.append(f"  - {f}: {val}")

    warnings = classification_result.get("warnings", [])
    if warnings:
        lines.append("\n Warnings:")
        for w in warnings:
            lines.append(f"  - {w}")

    tx = classification_result.get("treatment_suggestion")
    if tx:
        lines.append(f"\n Treatment suggestion: {tx}")

    if classification_result.get("hp_eradication_suggested"):
        lines.append("\n Suggestion: Helicobacter pylori eradication therapy is advised.(REC.38-40)")

    if data.get("smoking_status") == "smoker":
        lines.append("\n Suggestion: smokingcessation in individuals with precancerousconditions or after endoscopic treatment of superficial lesions. (REC.42).")

    missing = classification_result.get("missing_fields", [])
    if missing:
        lines.append(f"\n Missing fields (may affect decision completeness): {', '.join(missing)}")

    return "\n".join(lines)