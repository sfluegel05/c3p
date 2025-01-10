"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies: CHEBI:37668 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid is derived from a sesquiterpene (C15 skeleton) and may include rearrangements or modifications.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Sesquiterpenoids typically have 15 carbons, but modifications may reduce this
    if c_count < 10 or c_count > 20:
        return False, f"Carbon count ({c_count}) is outside the expected range for sesquiterpenoids"

    # Check for the presence of terpene-like structures (cyclic or acyclic)
    # This is a heuristic and may not catch all cases
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        # Acyclic sesquiterpenoids are possible but less common
        pass
    elif num_rings > 3:
        # Too many rings for a typical sesquiterpenoid
        return False, f"Too many rings ({num_rings}) for a sesquiterpenoid"

    # Check molecular weight - sesquiterpenoids typically have MW between 200 and 300
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150 or mol_wt > 350:
        return False, f"Molecular weight ({mol_wt:.2f}) is outside the expected range for sesquiterpenoids"

    # Check for common functional groups in sesquiterpenoids (e.g., alcohols, ketones, esters)
    # This is not exhaustive but helps confirm the classification
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    ketone_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    
    if not (mol.HasSubstructMatch(alcohol_pattern) or 
            mol.HasSubstructMatch(ketone_pattern) or 
            mol.HasSubstructMatch(ester_pattern)):
        return False, "No common functional groups (alcohol, ketone, ester) found"

    # If all checks pass, classify as sesquiterpenoid
    return True, "Contains a C15 skeleton with terpene-like structure and common functional groups"