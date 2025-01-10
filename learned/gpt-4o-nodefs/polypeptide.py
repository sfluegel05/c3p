"""
Classifies: CHEBI:15841 polypeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.
    Polypeptides are long chains of amino acids with peptide bonds (-CO-NH-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for improved peptide bond pattern (-CO-NH-)
    peptide_bond_pattern = Chem.MolFromSmarts("N[C;R0][C;R0](=O)")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) < 2:
        return False, f"Detected {len(peptide_bond_matches)} peptide bonds, need at least 2 for a polypeptide"
    
    # Look for enhanced characteristic polypeptide backbone
    backbone_pattern = Chem.MolFromSmarts("N[C;!R]=O")
    backbone_matches = mol.GetSubstructMatches(backbone_pattern)
    if len(backbone_matches) < 3:
        return False, "No characteristic polypeptide backbone found"

    # Evaluate molecular complexity and weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    
    if mol_wt < 800 or num_rotatable_bonds < 5 or num_rings > 0:
        return False, "Molecular complexity too low for a polypeptide"

    return True, "Contains peptide bonds and polypeptide backbone structure with sufficient molecular complexity"