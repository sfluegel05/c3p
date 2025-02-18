"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: CHEBI:16670 polypeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide (10+ amino acid residues) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count all amide bonds using RDKit's built-in function
    amide_bonds = rdMolDescriptors.CalcNumAmideBonds(mol)
    
    # Check molecular weight to filter small molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # For linear peptides: residues = amide_bonds + 1
    # For cyclic peptides: residues = amide_bonds
    # Require at least 9 amide bonds (10 residues linear or 9 cyclic)
    # Also check molecular weight > 500 Da as heuristic
    if amide_bonds >= 9 and mol_wt >= 500:
        return True, f"Contains {amide_bonds} amide bonds, indicating 10 or more residues"
    else:
        return False, f"Only {amide_bonds} amide bonds and/or molecular weight too low ({mol_wt:.1f} Da)"