"""
Classifies: CHEBI:46895 lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide contains a peptide (amino acid chain) and a lipid (long hydrocarbon chain).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Peptide backbone: Look for repeating N-C-C(=O) pattern
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX4][CX4](=[OX1])")
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    if len(peptide_matches) < 2:
       return False, "Too few peptide backbone units found"

    # 2. Lipid component: Look for a long aliphatic chain (at least 4 carbons) linked to a carbonyl
    lipid_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(lipid_chain_pattern)
    
    if not chain_matches:
        return False, "No long aliphatic chain found"

    # 3. Check for linkage between lipid and peptide via a carbonyl
    linked_pattern_N = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX3](=[OX1])[NX3]")
    linked_pattern_O = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX3](=[OX1])[OX2][#6]")
    
    linked_matches_N = mol.GetSubstructMatches(linked_pattern_N)
    linked_matches_O = mol.GetSubstructMatches(linked_pattern_O)


    if not linked_matches_N and not linked_matches_O:
        return False, "No connectivity between peptide and lipid"
    
    # 4. Molecular weight check. Lipopeptides are larger molecules, typically.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for a lipopeptide"

    return True, "Contains a peptide chain linked to a lipid chain"