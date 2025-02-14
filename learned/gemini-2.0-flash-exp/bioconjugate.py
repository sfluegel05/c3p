"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate is defined as a molecular entity consisting of at least 2 biological molecules covalently linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check molecule size (molecular weight)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
         return False, "Molecular weight too low to be a bioconjugate"
    
    # 2. Check for diverse heteroatoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    
    if n_count + o_count + p_count + s_count < 3:
      return False, "Too few heteroatoms for bioconjugate"
    
    # 3. Search for common linkages
    amide_pattern = Chem.MolFromSmarts("N[CX3](=O)")
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    thioether_pattern = Chem.MolFromSmarts("C-S-C")
    disulfide_pattern = Chem.MolFromSmarts("S-S")
    
    linkage_matches = mol.GetSubstructMatches(amide_pattern)
    linkage_matches.extend(mol.GetSubstructMatches(ester_pattern))
    linkage_matches.extend(mol.GetSubstructMatches(thioether_pattern))
    linkage_matches.extend(mol.GetSubstructMatches(disulfide_pattern))

    if len(linkage_matches) < 1:
        return False, "No common linkages found"

    # 4. Check number of rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Too few rotatable bonds"

    # 5. check for common biomolecule patterns (simplified version)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=[OX1])[OX2]")
    nucleotide_pattern = Chem.MolFromSmarts("P(=O)(O)O[CX4]")

    if mol.HasSubstructMatch(amino_acid_pattern) or mol.HasSubstructMatch(nucleotide_pattern) or mol_wt > 500:
      return True, "Contains biomolecule fragments and appropriate linkages"

    return False, "Does not meet criteria of bioconjugate"