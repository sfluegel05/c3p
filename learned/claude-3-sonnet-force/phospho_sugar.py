"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: CHEBI:83830 phospho sugar
'Any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.'
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify monosaccharide backbone
    # Use a recursive SMARTS pattern to match common monosaccharide ring structures
    monosaccharide_pattern = Chem.MolFromSmarts("[OR2][CR2][CR1][OR2][CR2][CR1]")
    if not mol.HasSubstructMatch(monosaccharide_pattern):
        return False, "No monosaccharide backbone found"

    # Check for presence of phosphate group(s)
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # Check if phosphate group(s) are connected to the monosaccharide backbone
    for phosphate_match in phosphate_matches:
        phosphate_atom = mol.GetAtomWithIdx(phosphate_match[0])
        connected_atoms = [mol.GetAtomWithIdx(neighbor).GetSymbol() for neighbor in phosphate_atom.GetNeighbors()]
        if "O" in connected_atoms:
            # Phosphate group is connected to an oxygen atom, which could be part of the monosaccharide
            pass
        else:
            return False, "Phosphate group not connected to monosaccharide backbone"

    # Check for common functional groups or modifications
    amino_pattern = Chem.MolFromSmarts("N")
    has_amino_group = mol.HasSubstructMatch(amino_pattern)

    nucleobase_pattern = Chem.MolFromSmarts("c1nc[nH]c1")
    has_nucleobase = mol.HasSubstructMatch(nucleobase_pattern)

    # Additional checks or patterns for other common modifications can be added here

    if has_amino_group:
        return True, "Contains monosaccharide backbone with phosphate group(s) and amino group"
    elif has_nucleobase:
        return True, "Contains monosaccharide backbone with phosphate group(s) and nucleobase"
    else:
        return True, "Contains monosaccharide backbone with phosphate group(s)"