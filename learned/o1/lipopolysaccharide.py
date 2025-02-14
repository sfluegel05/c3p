"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: lipopolysaccharide
"""

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    A lipopolysaccharide consists of a lipid A moiety (N-acylated glucosamine disaccharide with fatty acids),
    a core polysaccharide (including KDO and heptoses), and an O-antigen polysaccharide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Initialize feature counters
    features_found = []
    criteria_needed = 3  # Number of criteria to meet
    criteria_met = 0

    # 1. Check for sugar rings (pyranoses and furanoses)
    sugar_smarts = '[OX2H][CX4H][CX4H][OX2][CX4H][CX4H][OX2H]'  # Simplified pattern for sugars
    sugar = Chem.MolFromSmarts(sugar_smarts)
    sugar_matches = mol.GetSubstructMatches(sugar)
    num_sugars = len(sugar_matches)

    if num_sugars >= 3:
        criteria_met += 1
        features_found.append(f"Contains {num_sugars} sugar units")
    else:
        features_found.append(f"Contains {num_sugars} sugar units, less than required")

    # 2. Check for long-chain fatty acids (C12 or longer)
    fatty_acid_smarts = 'C(=O)[O][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]'  # Ester linkage to C12 chain
    fatty_acid = Chem.MolFromSmarts(fatty_acid_smarts)
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid)
    num_fatty_acids = len(fatty_acid_matches)

    if num_fatty_acids >= 2:
        criteria_met +=1
        features_found.append(f"Contains {num_fatty_acids} long-chain fatty acid esters")
    else:
        features_found.append(f"Contains {num_fatty_acids} long-chain fatty acid esters, less than required")

    # 3. Check for glucosamine units
    glucosamine_smarts = '[C@@H]([CX4H][OX2H])[CX4H][OX2H][CX4H][OX2H][CX4H][OX2H][CX4H][NX3H2]'
    glucosamine = Chem.MolFromSmarts(glucosamine_smarts)
    glucosamine_matches = mol.GetSubstructMatches(glucosamine)
    num_glucosamines = len(glucosamine_matches)

    if num_glucosamines >= 1:
        criteria_met +=1
        features_found.append("Contains glucosamine units")
    else:
        features_found.append("No glucosamine units found")

    # 4. Check for KDO units (octulosonic acid)
    kdo_smarts = '[CX3](=O)[CX4H][OX2H][CX4H][OX2H][CX4H][OX2H][CX4H][OX2H][CX4][OX2][CX3](=O)[O-]'
    kdo = Chem.MolFromSmarts(kdo_smarts)
    kdo_matches = mol.GetSubstructMatches(kdo)
    num_kdo = len(kdo_matches)

    if num_kdo >= 1:
        criteria_met +=1
        features_found.append("Contains KDO (octulosonic acid) units")
    else:
        features_found.append("No KDO (octulosonic acid) units found")

    # 5. Check for heptose units
    heptose_smarts = '[C@H]1([OX2H])[C@H]([OX2H])[C@H]([OX2H])[C@H]([OX2H])[C@H]([OX2H])[C@H](O)[C@H]1[OX2H]'
    heptose = Chem.MolFromSmarts(heptose_smarts)
    heptose_matches = mol.GetSubstructMatches(heptose)
    num_heptose = len(heptose_matches)

    if num_heptose >= 2:
        criteria_met +=1
        features_found.append(f"Contains {num_heptose} heptose units")
    else:
        features_found.append(f"Contains {num_heptose} heptose units, less than required")

    # 6. Check for multiple ester and amide bonds
    num_esters = len(mol.GetSubstructMatches(Chem.MolFromSmarts('C(=O)O')))
    num_amides = len(mol.GetSubstructMatches(Chem.MolFromSmarts('C(=O)N')))
    total_bonds = num_esters + num_amides

    if total_bonds >= 5:
        criteria_met +=1
        features_found.append(f"Contains {total_bonds} ester and amide bonds")
    else:
        features_found.append(f"Contains {total_bonds} ester and amide bonds, less than required")

    # 7. Check Molecular Weight
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt >= 1000:
        criteria_met +=1
        features_found.append(f"Molecular weight is {mol_wt:.2f} Da")
    else:
        features_found.append(f"Molecular weight is {mol_wt:.2f} Da, less than required")

    # Determine if criteria are met
    if criteria_met >= criteria_needed:
        return True, "; ".join(features_found)
    else:
        return False, "; ".join(features_found)