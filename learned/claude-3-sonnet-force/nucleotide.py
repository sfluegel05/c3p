"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:36973 nucleotide
A nucleotide is a nucleoside phosphate resulting from the condensation of the 3 or 5 hydroxy group of a nucleoside with phosphoric acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for nucleoside base
    base_pattern = Chem.MolFromSmarts("[nr5]1[nr5]c2nc[nr5]c2[nr5]1")
    if not mol.HasSubstructMatch(base_pattern):
        return False, "No nucleoside base found"

    # Look for sugar moiety
    sugar_pattern = Chem.MolFromSmarts("[OX2r5]1[Cr5][Cr5]([Or5])[Cr5]([Or5])[Cr5]1")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moiety found"

    # Look for phosphate group(s)
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)([O-,O])")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # Check for connectivity between components
    base_atom = mol.GetSubstructMatch(base_pattern)[0]
    sugar_atom = mol.GetSubstructMatch(sugar_pattern)[0]
    phosphate_atom = phosphate_matches[0][0]

    # Check if base and sugar are connected
    base_sugar_bond = mol.GetBondBetweenAtoms(base_atom, sugar_atom)
    if base_sugar_bond is None:
        return False, "Nucleoside base and sugar moiety not connected"

    # Check if sugar and phosphate are connected
    sugar_phosphate_bond = mol.GetBondBetweenAtoms(sugar_atom, phosphate_atom)
    if sugar_phosphate_bond is None:
        return False, "Sugar moiety and phosphate group not connected"

    # Additional checks (optional)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 1000:
        return False, "Molecular weight outside typical range for nucleotides"

    return True, "Contains a nucleoside base, sugar moiety, and phosphate group(s) in the correct structural arrangement"