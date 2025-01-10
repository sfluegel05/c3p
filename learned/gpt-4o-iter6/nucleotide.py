"""
Classifies: CHEBI:36976 nucleotide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide consists of a nucleoside (a nitrogenous base linked to a sugar) bonded to a phosphate group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a nucleotide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify at least one phosphate group, considering variations (mono, di, tri, cyclic)
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX1])[OX2]")  # Adapt more generic pattern to identify phosphate groups
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_matches < 1:
        return False, "No phosphate group found"

    # Check for the presence of ribose or deoxyribose sugar rings
    ribose_pattern = Chem.MolFromSmarts("OC1COC(O)C(O)C1[OX2H]")  # Adjusted pattern to ensure correct stereochemistry
    deoxyribose_pattern = Chem.MolFromSmarts("OC1COC(O)C1[OX2H]")

    if not mol.HasSubstructMatch(ribose_pattern) and not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No compatible sugar ring (ribose or deoxyribose) found"

    # Check for nitrogenous bases or any common nucleotide bases
    base_patterns = [
        Chem.MolFromSmarts("n1cnc2c1ncnc2"),  # Adequate for purine ring systems
        Chem.MolFromSmarts("n1c[nH]c2c1ncnc2"),  # Adenine, Guanine
        Chem.MolFromSmarts("c1c[nH]c(=O)[nH]c1=O"),  # Uracil, Thymine
        Chem.MolFromSmarts("c1[nH]c(=O)[nH]c2c1ncnc2"),  # Cytosine variation
    ]
    
    base_match_found = any(mol.HasSubstructMatch(base_pattern) for base_pattern in base_patterns)
    if not base_match_found:
        return False, "No nitrogenous base found in the structure"

    return True, "Contains nucleoside with phosphate, sugar, and base attached appropriately"