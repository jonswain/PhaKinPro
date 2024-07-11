"""PhaKinPro main functionality"""

import _pickle as cPickle
import argparse
import bz2
import glob
import gzip
import io
import os
from io import StringIO
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rdkit
from config import AD_DICT, CLASSIFICATION_DICT, MODEL_DICT, MODEL_DICT_INVERT
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MolFromSmiles
from rdkit.Chem.Draw import SimilarityMaps
from scipy.spatial.distance import cdist
from sklearn.base import BaseEstimator
from tqdm import tqdm


def run_prediction(
    model: BaseEstimator,
    model_data: dict,
    fingerprints: pd.DataFrame,
    calculate_ad: bool = True,
) -> tuple[list, np.array, list | None]:
    """Makes prediction using a model and the fingerprints.

    Args:
        model (BaseEstimator): The model to use for prediction.
        model_data (dict): The data used to train the model.
        fingerprints (pd.DataFrame): Fingerprints to predict on
        calculate_ad (bool, optional): Calculate applicability domain. Defaults to True.

    Returns:
        tuple[list, np.array, list | None]: Predictions, probabilities, and inside AD.
    """
    pred_proba = model.predict_proba(fingerprints)[:, 1]
    preds = [1 if prob > 0.5 else 0 for prob in pred_proba]

    for idx, pred in enumerate(preds):
        if pred == 0:
            pred_proba[idx] = 1 - pred_proba[idx]

    if calculate_ad:
        ads = []
        for idx, row in fingerprints.iterrows():
            ad = model_data["D_cutoff"] > np.min(
                cdist(
                    model_data["Descriptors"].to_numpy(), np.array(row).reshape(1, -1)
                )
            )
            ads.append(ad)
        return preds, pred_proba, ads
    return preds, pred_proba, None


# TODO: Make this function accept all fingerprints and make predictions at once
def get_prob_map(model, smiles):
    def get_fp(mol, idx):
        fps = np.zeros((2048, 1))
        _fps = SimilarityMaps.GetMorganFingerprint(mol, idx, radius=3, nBits=2048)
        DataStructs.ConvertToNumpyArray(_fps, fps)
        return fps

    def get_proba(fps):
        return float(model.predict_proba(fps.reshape(1, -1))[:, 1])

    mol = Chem.MolFromSmiles(smiles)
    fig, _ = SimilarityMaps.GetSimilarityMapForModel(mol, get_fp, get_proba)
    imgdata = io.StringIO()
    fig.savefig(imgdata, format="svg")
    imgdata.seek(0)  # rewind the data
    plt.savefig(imgdata, format="svg", bbox_inches="tight")

    return imgdata.getvalue()


def multiclass_ranking(ordered_preds: list) -> int:
    idx = 0
    one_detected = False
    for i, o in enumerate(ordered_preds):
        if int(o) == 1:
            if not one_detected:
                idx = i + 1
                one_detected = True
        if int(o) == 0:
            if one_detected:
                idx = 0
                return idx
    return idx if idx != 0 else len(ordered_preds) + 1


def get_prediction_data(
    fingerprints: pd.DataFrame,
    calculate_ad: bool = True,
    make_prop_img: bool = False,
    **kwargs,
) -> list:

    # This checks if the found models are in the directory
    def default(key, d):
        if key in d.keys():
            return d[key]
        else:
            return False

    # Find all models
    models = sorted(
        [
            f
            for f in glob.glob(
                os.path.join(
                    os.path.dirname(os.path.realpath(__file__)), "./models/*.pgz"
                )
            )
        ],
        key=lambda x: x.split("_")[1],
    )
    # Find all data
    models_data = sorted(
        [
            f
            for f in glob.glob(
                os.path.join(
                    os.path.dirname(os.path.realpath(__file__)), "./models/*.pbz2"
                )
            )
        ],
        key=lambda x: x.split("_")[1],
    )

    values: dict = {}

    for model_endpoint, model_data in tqdm(
        zip(models, models_data), total=min(len(models), len(models_data))
    ):
        if not default(MODEL_DICT_INVERT[os.path.basename(model_endpoint)], MODEL_DICT):
            continue
        with gzip.open(model_endpoint, "rb") as f:
            model = cPickle.load(f)

        with bz2.BZ2File(model_data, "rb") as f:
            model_data = cPickle.load(f)

        preds, pred_probas, ads = run_prediction(
            model, model_data, fingerprints, calculate_ad=calculate_ad
        )

        if ads is None:
            ads = [None] * len(preds)

        # TODO: Make this work with fingerprints
        svg_str = ""
        # if make_prop_img:
        #     svg_str = get_prob_map(model, smiles

        values.setdefault(
            MODEL_DICT_INVERT[os.path.basename(model_endpoint)], []
        ).append(
            {
                "preds": preds,
                "probability": [
                    str(round(float(pred_proba) * 100, 2)) + "%"
                    for pred_proba in pred_probas
                ],
                "inside_ad": ads,
                "svg_str": svg_str,
            }
        )

    processed_results = []

    for idx, smiles in enumerate(fingerprints.index):
        processed_result = []
        for key, value in values.items():
            if key in [
                "Hepatic Stability",
                "Renal Clearance",
                "Plasma Half-life",
                "Oral Bioavailability",
            ]:
                val = []
                for model in value:
                    val.append(
                        [
                            model["preds"][idx],
                            model["probability"][idx],
                            model["inside_ad"][idx],
                            model["svg_str"],
                        ]
                    )

                new_pred = multiclass_ranking([_[0] for _ in val])
                if new_pred == 0:
                    processed_result.append(
                        [
                            key,
                            "Inconsistent result: no prediction",
                            "Very unconfident",
                            "NA",
                            "",
                        ]
                    )
                else:
                    # this is because of how the hierarchical model works
                    if new_pred in [1, 2]:
                        p = 0
                    else:
                        p = new_pred - 2
                    processed_result.append(
                        [
                            key,
                            CLASSIFICATION_DICT[key][new_pred],
                            val[p][1],
                            val[p][2],
                            val[p][3],
                        ]
                    )
            else:
                processed_result.append(
                    [
                        key,
                        CLASSIFICATION_DICT[key][val[0][0]],
                        val[0][1],
                        val[0][2],
                        val[0][3],
                    ]
                )

        processed_results.append((smiles, processed_result))

    return processed_results


def create_results_frame(fingerprints, calculate_ad=False) -> pd.DataFrame:
    # Create columns for results
    columns = []
    for _key in MODEL_DICT.keys():
        columns.append(_key)
        columns.append(_key + "_proba")
        if calculate_ad:
            columns.append(_key + "_AD")

    # Get results
    results = get_prediction_data(fingerprints, calculate_ad=calculate_ad, **MODEL_DICT)

    final_table = []
    for smiles, data in results:
        cols = {}
        for model_name, pred, pred_proba, ad, _ in data:
            row = {}
            try:
                pred_proba = float(pred_proba[:-1]) / 100  # covert back to 0-1 float
                row[model_name] = pred
                row[model_name + "_proba"] = (
                    pred_proba if pred == 1 else 1 - pred_proba
                )  # this is to make sure its proba for class 1
            except ValueError:
                row[model_name] = "No prediction"  # if pred_proba is string skip
            if calculate_ad:
                row[model_name + "_AD"] = ad

            cols.update(row)
        final_table.append({smiles: cols})

    rows = [pd.DataFrame.from_dict(row, orient="index") for row in final_table]
    return pd.concat(rows)


def create_fingerprints(mols_list: list[rdkit.Chem.rdchem.Mol | None]) -> pd.DataFrame:
    """Create Morgan fingerprints for a list of SMILES strings.

    Args:
        mols_list (list[rdkit.Chem.rdchem.Mol | None]): List of RDKit molecules.

    Returns:
        pd.DataFrame: DataFrame of fingerprints.
    """
    fps = []
    for mol in tqdm(mols_list):
        fp = np.zeros((2048, 1))
        _fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=2048)
        DataStructs.ConvertToNumpyArray(_fp, fp)
        fps.append(fp.astype(int))
    return pd.DataFrame(fps, columns=[f"Bit_{i}" for i in range(2048)])


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--infile", type=str, required=True, help="location to csv of SMILES"
    )
    parser.add_argument(
        "--outfile",
        type=str,
        default=os.path.join(os.getcwd(), "phakin_output.csv"),
        help="location and file name for output",
    )
    parser.add_argument(
        "--smiles_col",
        type=str,
        default="SMILES",
        help="column name containing SMILES of interest",
    ),
    parser.add_argument(
        "--id_col",
        type=str,
        help="column name containing compound identifiers",
    ),
    parser.add_argument("--ad", action="store_true", help="calculate the AD")
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    """Main function to run PhaKinPro predictions."""
    # Validate input file
    if not Path(args.infile).exists():
        raise FileNotFoundError(f"Input file {args.infile} not found.")
    input_df = pd.read_csv(Path(args.infile))
    if args.smiles_col not in input_df.columns:
        raise ValueError(f"Column {args.smiles_col} not found in {args.infile}")
    if args.id_col not in input_df.columns:
        raise ValueError(f"Column {args.smiles_col} not found in {args.infile}")

    # Get SMILES and identifiers from CSV input
    smiles = input_df[args.smiles_col]
    identifiers = input_df[args.id_col]

    # Create predictions dataframe
    prediction_df = pd.DataFrame({"smiles": smiles})
    if args.id_col:
        prediction_df[args.id_col] = identifiers
    print("Checking input SMILES...")
    prediction_df["mols"] = [MolFromSmiles(s) for s in smiles]
    invalid_smiles = prediction_df[prediction_df["mols"].isnull()]
    if len(invalid_smiles) > 0:
        print(f"Invalid SMILES at indices: {invalid_smiles.index.to_list()}")
        # Rename invalid SMILES f"(invalid){x}"
        prediction_df.loc[invalid_smiles.index, "smiles"] = (
            "(invalid) " + prediction_df.loc[invalid_smiles.index, "smiles"]
        )

    # Create fingerprints
    print("Creating fingerprints...")
    fingerprints = create_fingerprints(
        prediction_df[prediction_df["mols"].notnull()]["mols"]
    )
    fingerprints.index = prediction_df[prediction_df["mols"].notnull()]["smiles"]

    # Make predictions
    print("Making predictions...")
    predictions = create_results_frame(fingerprints)
    predictions.index.name = "smiles"

    # Merge predictions with input data, tidy, and save
    merged_df = prediction_df.merge(
        predictions.reset_index(), left_on="smiles", right_on="smiles", how="left"
    )

    merged_df.rename(
        columns={"smiles": args.smiles_col, "identifier": args.id_col}, inplace=True
    )
    merged_df.drop("mols", axis=1).to_csv(Path(args.outfile), index=False)
    print("Finished! Results saved to:", args.outfile)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
