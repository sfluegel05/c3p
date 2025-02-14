"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: CHEBI:35930 octanoate ester
Any fatty acid ester in which the carboxylic acid component is octanoic acid (caprylic acid).
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for octanoic acid (caprylic acid) substructure
    octanoic_acid_pattern = Chem.MolFromSmarts("CCCCCCCC(=O)O")
    octanoic_acid_matches = mol.GetSubstructMatches(octanoic_acid_pattern)
    
    # Look for any ester bond (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Check molecular weight and atom counts
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if len(octanoic_acid_matches) > 0 and len(ester_matches) > 0 and mol_wt > 200 and c_count >= 8 and o_count >= 2:
        return True, "Contains octanoic acid component and ester bond"
    else:
        return False, "Does not meet the criteria for an octanoate ester"