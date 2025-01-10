"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are characterized by a heterocyclic ring (often benzene) with oxygen-containing functional groups,
    and often have long hydrocarbon chains or amide/ester linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a benzene ring (common in cannabinoids)
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene_pattern):
        return False, "No benzene ring found"

    # Check for oxygen-containing functional groups (e.g., hydroxyl, ester, amide)
    oxygen_pattern = Chem.MolFromSmarts("[OX2]")
    oxygen_matches = mol.GetSubstructMatches(oxygen_pattern)
    if len(oxygen_matches) < 1:
        return False, "No oxygen-containing functional groups found"

    # Check for long hydrocarbon chains or amide/ester linkages
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    
    if not (mol.HasSubstructMatch(long_chain_pattern) or 
            mol.HasSubstructMatch(amide_pattern) or 
            mol.HasSubstructMatch(ester_pattern)):
        return False, "No long hydrocarbon chains or amide/ester linkages found"

    # Check molecular weight - cannabinoids typically have a higher molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for a cannabinoid"

    return True, "Contains a benzene ring with oxygen-containing functional groups and long hydrocarbon chains or amide/ester linkages"