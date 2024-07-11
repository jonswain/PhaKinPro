"""Config file for PhaKinPro"""

MODEL_DICT = {
    "Hepatic Stability": [
        "Dataset_01B_hepatic-stability_15min_imbalanced-morgan_RF.pgz",
        "Dataset_01C_hepatic-stability_30min_imbalanced-morgan_RF.pgz",
        "Dataset_01D_hepatic-stability_60min_imbalanced-morgan_RF.pgz",
    ],
    "Microsomal Half-life Sub-cellular": [
        "Dataset_02A_microsomal-half-life-subcellular_imbalanced-morgan_RF.pgz"
    ],
    "Microsomal Half-life Tissue": [
        "Dataset_02B_microsomal-half-life_30-min_binary_unbalanced_morgan_RF.pgz"
    ],
    "Renal Clearance": [
        "dataset_03_renal-clearance_0.1-threshold_balanced-morgan_RF.pgz",
        "dataset_03_renal-clearance_0.5-threshold_imbalanced-morgan_RF.pgz",
        "dataset_03_renal-clearance_1.0-threshold_balanced-morgan_RF.pgz",
    ],
    "BBB Permeability": ["dataset_04_bbb-permeability_balanced-morgan_RF.pgz"],
    "CNS Activity": ["dataset_04_cns-activity_1464-compounds_imbalanced-morgan_RF.pgz"],
    "CACO2": ["Dataset_05A_CACO2_binary_unbalanced_morgan_RF.pgz"],
    "Plasma Protein Binding": [
        "Dataset_06_plasma-protein-binding_binary_unbalanced_morgan_RF.pgz"
    ],
    "Plasma Half-life": [
        "Dataset_08_plasma_half_life_12_hr_balanced-morgan_RF.pgz",
        "Dataset_08_plasma_half_life_1_hr_balanced-morgan_RF.pgz",
        "Dataset_08_plasma_half_life_6_hr_imbalanced-morgan_RF.pgz",
    ],
    "Microsomal Intrinsic Clearance": [
        "Dataset_09_microsomal-intrinsic-clearance_12uL-min-mg-threshold-imbalanced-morgan_RF.pgz"
    ],
    "Oral Bioavailability": [
        "dataset_10_oral_bioavailability_0.5_threshold_imbalanced-morgan_RF.pgz",
        "dataset_10_oral_bioavailability_0.8_balanced-morgan_RF.pgz",
    ],
}

# This is to seaerch the keys of the dict for the model paths and return the model name.
MODEL_DICT_INVERT = {v: key for key, val in MODEL_DICT.items() for v in val}

CLASSIFICATION_DICT = {
    "Hepatic Stability": {
        1: "Hepatic stability <= 50% at 15 minutes",
        2: "Hepatic stability <= 50% between 15 and 30 minutes",
        3: "Hepatic stability <= 50% between 30 and 60 minutes",
        4: "Hepatic stability > 50% at 60 minutes",
    },
    "Microsomal Half-life Sub-cellular": {
        0: "Sub-cellular Hepatic Half-life > 30 minutes",
        1: "Sub-cellular Hepatic Half-life <= 30 minutes",
    },
    "Microsomal Half-life Tissue": {
        0: "Tissue Hepatic Half-life > 30 minutes",
        1: "Tissue Hepatic Half-life <= 30 minutes",
    },
    "Renal Clearance": {
        1: "Renal clearance below 0.10 ml/min/kg",
        2: "Renal clearance between 0.10 and 0.50 ml/min/kg",
        3: "Renal clearance between 0.50 and 1.00 ml/min/kg",
        4: "Renal clearance above 1.00 ml/min/kg",
    },
    "BBB Permeability": {
        0: "Does not permeate blood brain barrier",
        1: "Does permeate blood brain barrier",
    },
    "CNS Activity": {
        0: "Does not exhibit central nervous system activity",
        1: "Does exhibit central nervous system activity",
    },
    "CACO2": {0: "Does not permeate Caco-2", 1: "Does permeate Caco-2"},
    "Plasma Protein Binding": {
        0: "Plasma protein binder",
        1: "Weak/non plasma protein binder",
    },
    "Plasma Half-life": {
        1: "Half-life below 1 hour",
        2: "Half-life between 1 and 6 hours",
        3: "Half-life between 6 and 12 hours",
        4: "Half-life above 12 hours",
    },
    "Microsomal Intrinsic Clearance": {
        0: "Microsomal intrinsic clearance < 12 uL/min/mg",
        1: "Microsomal intrinsic clearance >= 12 uL/min/mg",
    },
    "Oral Bioavailability": {
        1: "Less than 0.5 F",
        2: "Between 0.5 and 0.8 F",
        3: "Above 0.8 F",
    },
}

AD_DICT = {True: "Inside", False: "Outside"}
