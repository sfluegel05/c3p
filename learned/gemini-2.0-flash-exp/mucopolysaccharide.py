"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
import numpy as np

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is likely a mucopolysaccharide based on its SMILES string using heuristics and a scoring system.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) - True if likely a mucopolysaccharide, False otherwise, along with a reason.
               Returns (None, None) if SMILES processing fails.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None

    score = 0
    reasons = []

    # Define SMARTS patterns (more specific)
    # Uronic acid - modified to accept 2- or 3- modified uronic acids
    uronic_acid_pattern = Chem.MolFromSmarts("[C;R1][C;R1]([O])[C;R1]([O])[C;R1][C;R1][C](=O)[O;H0]")
    uronic_acid_matches = mol.GetSubstructMatches(uronic_acid_pattern)

    # Sulfate group (more specific)
    sulfate_pattern = Chem.MolFromSmarts("O[S](=[O])(=[O])[O]")
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)

    # Amino group in a sugar ring (more specific)
    amino_pattern = Chem.MolFromSmarts("[NX3;R1]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)

    # Glycosidic linkages (O-C-O within a ring)
    glycosidic_pattern = Chem.MolFromSmarts("[O;R1][C;R1][O;R1]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)

    # Basic ring oxygens
    oxygen_pattern = Chem.MolFromSmarts("[OX2;R1]")
    oxygen_matches = mol.GetSubstructMatches(oxygen_pattern)

    # Count carbons and non-standard atoms
    carbon_pattern = Chem.MolFromSmarts("[CX4]")
    carbon_matches = mol.GetSubstructMatches(carbon_pattern)

    other_atoms_count = mol.GetNumAtoms() - (sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
                                            + sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
                                            + sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
                                            + sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
                                            + sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16))

    n_atoms = mol.GetNumAtoms()
    if n_atoms < 20:
        return False, "Too small for mucopolysaccharide"

    # Score based on detected features
    if len(uronic_acid_matches) > 0:
        score += 2
        reasons.append(f"Found {len(uronic_acid_matches)} uronic acid units")
    if len(sulfate_matches) > 0:
        score += 2
        reasons.append(f"Found {len(sulfate_matches)} sulfate groups")
    if len(amino_matches) > 0:
      score += 1
      reasons.append(f"Found {len(amino_matches)} amino groups")
    if len(glycosidic_matches) > 2:
        score += 1
        reasons.append(f"Found {len(glycosidic_matches)} glycosidic linkages")

    #Size-related conditions (more strict)
    if len(oxygen_matches) < 4:
      return False, "Too few oxygen atoms"
    if len(carbon_matches) < 10:
        return False, "Too few carbon atoms"
    if len(oxygen_matches) > 2 * len(carbon_matches):
        return False, "Too many oxygen atoms"
    if other_atoms_count > 5:
        return False, "Too many non-standard atoms"

    # Molecular weight and rotatable bonds
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a mucopolysaccharide"

    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5 :
      return False, "Too few rotatable bonds for a mucopolysaccharide"


    if score < 2 :
      return False, f"Low score {score}. {';'.join(reasons)}."


    return True, f"Likely mucopolysaccharide: score {score}. {';'.join(reasons)}"