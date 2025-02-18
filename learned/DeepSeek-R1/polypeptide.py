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
    
    # SMARTS pattern for peptide bonds in backbone: CONH connected to alpha carbon
    peptide_bond = Chem.MolFromSmarts('[CX3](=O)[NX3H1][CX4H2]')
    matches = mol.GetSubstructMatches(peptide_bond)
    backbone_amide_count = len(matches)
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Check criteria: at least 9 backbone amide bonds (linear: 10 residues) or 10 (cyclic: 10 residues)
    # and molecular weight over 500 as heuristic filter
    if backbone_amide_count >= 9 and mol_wt >= 500:
        return True, f"Contains {backbone_amide_count} backbone amide bonds, indicating 10 or more residues"
    else:
        return False, (f"Only {backbone_amide_count} backbone amide bonds and/or molecular weight too low "
                      f"({mol_wt:.1f} Da)")