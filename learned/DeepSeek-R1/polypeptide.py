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
    
    # Define peptide bond SMARTS: C(=O)-N connected to a carbon
    peptide_smarts = Chem.MolFromSmarts('[CX3](=[OX1])-[NX3][CX4]')
    peptide_bonds = len(mol.GetSubstructMatches(peptide_smarts))
    
    # For linear peptides: residues = peptide_bonds + 1
    # For cyclic peptides: residues = peptide_bonds
    # Check if either case gives >=10 residues
    if peptide_bonds >= 9:  # Linear (9+1=10) or cyclic (>=10)
        return True, f"Contains {peptide_bonds} peptide bonds, indicating 10 or more residues"
    else:
        return False, f"Only {peptide_bonds} peptide bonds, insufficient for 10 residues"