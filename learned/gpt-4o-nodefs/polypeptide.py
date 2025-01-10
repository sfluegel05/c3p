"""
Classifies: CHEBI:15841 polypeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.
    Polypeptides are characterized by sequences of peptide bonds (-CO-NH-) and amino acid structures.

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

    # The peptide bond pattern may need refinement to capture all necessary forms, including cyclized forms
    peptide_bond_pattern = Chem.MolFromSmarts("N[C;R0][C;R0](=O)")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) < 5:  # Assuming at least 5 peptide bonds indicate a polypeptide
        return False, f"Detected {len(peptide_bond_matches)} peptide bonds, need at least 5 for a polypeptide"

    # Count the number of nitrogen atoms generally found in polypeptides
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if nitrogen_count < 10:  # Expecting at least 10 nitrogen atoms typically found in short to long peptides
        return False, f"Detected {nitrogen_count} nitrogen atoms, too low for a polypeptide"

    # Check molecular weight to ensure we capture larger polypeptide molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:  # Molecular weight of 500 Dalton is typically a minimum for a polypeptide
        return False, f"Molecular weight too low ({mol_wt} Da) for a polypeptide"

    # Make sure there's a continuous backbone structure, i.e., no isolated segments of peptide-like groups
    molecular_fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(molecular_fragments) > 1:
        return False, "Polypeptide structure should not be fragmented"

    return True, "Contains characteristic peptide bonds and nitrogen-rich backbone indicative of a polypeptide"