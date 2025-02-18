"""
Classifies: CHEBI:33839 macromolecule
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    Considers molecular weight, repeating structural motifs, and complexity.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macromolecule, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight of {mol_wt:.2f} Da is below typical macromolecule threshold."

    # Identify repeating structural motifs
    # Check for peptide bonds, glycosidic linkages, nuclear acid bases
    peptide_bond = Chem.MolFromSmarts('C(=O)N')
    glycosidic_bond = Chem.MolFromSmarts('O[C@@H]1C(O)C1')
    nucleic_acid_base = Chem.MolFromSmarts('C1=NC=CN1')

    # Match counts
    peptide_count = len(mol.GetSubstructMatches(peptide_bond))
    glycosidic_count = len(mol.GetSubstructMatches(glycosidic_bond))
    nucleic_acid_count = len(mol.GetSubstructMatches(nucleic_acid_base))

    # Complex and multi-pattern structures correspond to macromolecules
    if peptide_count >= 5 or glycosidic_count >= 3 or nucleic_acid_count >= 2:
        return True, "Recurring polymeric structures suggest macromolecular nature."
    
    # Additional complexity checks
    if mol.GetNumAtoms() > 100:
        return True, "Large molecule, suggesting macromolecular structure."

    return False, "Could not identify distinct macromolecular features."