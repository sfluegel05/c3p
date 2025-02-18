"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: CHEBI:33839 macromolecule
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule is a high molecular weight molecule with repeating structural units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macromolecule, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular weight
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low"

    # Check for protein-like structure (multiple amide bonds)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    amide_count = len(mol.GetSubstructMatches(amide_pattern))
    if amide_count >= 5:
        return True, f"Contains {amide_count} amide bonds, indicative of a protein"

    # Check for polysaccharide-like structure (multiple ring ether linkages)
    # This SMARTS approximates ether oxygen connected to two carbons in rings
    glycosidic_pattern = Chem.MolFromSmarts("[C&r][O&!$(*-C=O)][C&r]")
    glycosidic_count = len(mol.GetSubstructMatches(glycosidic_pattern))
    if glycosidic_count >= 4:
        return True, f"Contains {glycosidic_count} glycosidic-like linkages"

    # Check for very high molecular weight as fallback
    if mol_wt > 5000:
        return True, f"Very high molecular weight ({mol_wt:.1f} Da)"

    return False, "Does not meet macromolecule criteria"