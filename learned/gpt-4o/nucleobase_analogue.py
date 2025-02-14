"""
Classifies: CHEBI:67142 nucleobase analogue
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is defined as a molecule that can substitute for a normal nucleobase in nucleic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for nitrogen-containing heterocycles
    pyrimidine_pattern = Chem.MolFromSmarts("c1ncnc1")  # A pyrimidine-like structure
    purine_pattern = Chem.MolFromSmarts("c1ncnc2[nH]ncc12")  # A purine-like structure

    if mol.HasSubstructMatch(pyrimidine_pattern) or mol.HasSubstructMatch(purine_pattern):
        return True, "Contains nitrogen-containing heterocycle similar to nucleobases"

    # Look for key functional groups common to nucleobases
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")  # Carbonyl group
    amino_pattern = Chem.MolFromSmarts("[NX3][H]")  # Amino group

    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    amino_matches = mol.GetSubstructMatches(amino_pattern)

    if len(carbonyl_matches) > 0 or len(amino_matches) > 0:
        return True, "Contains functional groups typical of nucleobase analogues"

    # Additional checks for modified nucleobases
    n_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    hetero_atm_count = rdMolDescriptors.CalcNumHeteroatoms(mol)

    if n_atoms > 1 and hetero_atm_count >= 2:
        return True, "Structure suggests a modified nucleobase"

    return False, "Does not match features of known nucleobase analogues"