"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: CHEBI:33697 macromolecule
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule is a molecule of high relative molecular mass, the structure of which
    essentially comprises the multiple repetition of units derived, actually or conceptually,
    from molecules of low relative molecular mass.

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
    
    # Check molecular weight - macromolecules typically > 1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, "Molecular weight too low for macromolecule"
    
    # Look for repeating substructures (peptides, polysaccharides, polyketides, etc.)
    peptide_pattern = Chem.MolFromSmarts("[N;X3]-[C;X3](-[C;X3])=[O;X1]")
    saccharide_pattern = Chem.MolFromSmarts("[O;X2]C[C;X4][O;X2]")
    polyketide_pattern = Chem.MolFromSmarts("[C;X3](=[O;X1])[C;X3](-[C;X3])=[C;X3]")
    
    if not (mol.HasSubstructMatch(peptide_pattern) or mol.HasSubstructMatch(saccharide_pattern) or mol.HasSubstructMatch(polyketide_pattern)):
        return False, "No repeating substructures found"
    
    # Check for presence of specific functional groups or structural motifs
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 3:
        return False, "Insufficient ring structures for macromolecule"
    
    hetero_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 1)
    if hetero_count < 5:
        return False, "Insufficient heteroatoms for macromolecule"
    
    # If all checks pass, classify as macromolecule
    return True, "Molecule has high molecular weight, repeating substructures, and structural features typical of macromolecules"