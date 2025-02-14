"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:36973 nucleotide
A nucleotide is a nucleoside phosphate resulting from the condensation of the 3 or 5 hydroxy group of a nucleoside with phosphoric acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem

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
    bases = ['Adenine', 'Guanine', 'Cytosine', 'Thymine', 'Uracil']
    base_found = False
    for base in bases:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(f'n{base}')):
            base_found = True
            break

    # Look for alternative/modified nucleoside bases
    if not base_found:
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic() and atom.GetSymbol() == 'N':
                env = Chem.FindAtomEnvironmentOfRadiusN(mol, 3, atom.GetIdx())
                if env.getIsNucleicAcidBase():
                    base_found = True
                    break

    if not base_found:
        return False, "No nucleoside base found"

    # Look for sugar moiety
    sugar_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetIsAromatic():
            env = Chem.FindAtomEnvironmentOfRadiusN(mol, 3, atom.GetIdx())
            if env.getIsSugar():
                sugar_found = True
                break

    if not sugar_found:
        return False, "No sugar moiety found"

    # Look for phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)([O-,O])")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # Additional checks for other structural features (if needed)
    # ...

    return True, "Contains a nucleoside base, sugar moiety, and phosphate group(s)"