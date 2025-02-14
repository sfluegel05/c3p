"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: CHEBI:37563 nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate is a ribosyl or deoxyribosyl derivative of a pyrimidine or purine base
    in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for pyrimidine or purine base
    base_pattern = Chem.MolFromSmarts("[*]1[*]c2[nH]c[nH]c2[*]c1[*]") # Purine
    base_pattern = base_pattern.Union(Chem.MolFromSmarts("[*]1[*]c[nH]c([nH]c1[*])[*]")) # Pyrimidine
    if not mol.HasSubstructMatch(base_pattern):
        return False, "No purine or pyrimidine base found"
    
    # Look for ribose/deoxyribose ring
    ribose_pattern = Chem.MolFromSmarts("[OX2]C1C(C(C(O[CX4H2]([OX2H0P](=O)(O)O)O1)O)O)O") # Ribose with phosphate
    deoxyribose_pattern = Chem.MolFromSmarts("[OX2]C1C(C(C(O[CX4H2]([OX2H0P](=O)(O)O)O1)O)N)") # Deoxyribose with phosphate
    if not (mol.HasSubstructMatch(ribose_pattern) or mol.HasSubstructMatch(deoxyribose_pattern)):
        return False, "No ribose or deoxyribose ring with phosphate found"
    
    # Check for additional phosphate groups (di-, tri-, tetra-phosphorylation)
    phosphate_pattern = Chem.MolFromSmarts("[OX2H0P](=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) > 4:
        return False, "More than 4 phosphate groups found"
    
    # Check molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1000:
        return False, "Molecular weight outside expected range for nucleoside 5'-phosphate"
    
    # Check atom counts
    n_atoms = mol.GetNumAtoms()
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if n_atoms < 20 or c_count < 10 or n_count < 5 or o_count < 5 or p_count < 1:
        return False, "Atom counts outside expected ranges for nucleoside 5'-phosphate"
    
    return True, "Molecule contains a purine or pyrimidine base attached to a ribose or deoxyribose ring with mono-, di-, tri- or tetra-phosphorylation at C-5"