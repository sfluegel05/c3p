"""
Classifies: CHEBI:17334 penicillin
"""
"""
Classifies: CHEBI:35489 penicillin
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    A penicillin is defined as any member of the group of substituted penams 
    containing two methyl substituents at position 2, a carboxylate substituent
    at position 3 and a carboxamido group at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for penam ring system
    penam_pattern = Chem.MolFromSmarts("[$(C2(N)C(=O)N(C1([CH3])([CH3])SC1(C)C)C2=O)]")
    if not mol.HasSubstructMatch(penam_pattern):
        return False, "No penam ring system found"
    
    # Check for methyl groups at position 2
    methyl_pattern = Chem.MolFromSmarts("[C@@](N)(C(=O)N1C(=O)[C@]2([CH3])([CH3])SC[C@@H]12)(C)")
    if not mol.HasSubstructMatch(methyl_pattern):
        return False, "Missing methyl groups at position 2"
    
    # Check for carboxylate group at position 3
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)([O-])")
    carboxylate_match = mol.GetSubstructMatches(carboxylate_pattern)
    if len(carboxylate_match) != 1:
        return False, "Missing or incorrect carboxylate group at position 3"
    
    # Check for carboxamido group at position 6
    carboxamido_pattern = Chem.MolFromSmarts("C(=O)N")
    carboxamido_match = mol.GetSubstructMatches(carboxamido_pattern)
    if len(carboxamido_match) != 1:
        return False, "Missing or incorrect carboxamido group at position 6"
    
    # Check for typical molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 600:
        return False, "Molecular weight outside typical range for penicillins"
    
    return True, "Contains penam ring system with two methyl groups at position 2, carboxylate at 3, and carboxamido at 6"

# Additional examples of true positives
print(is_penicillin("CC1(C)SC2N(C(=O)C(NC(=O)C(N)C3=CC=CC=C3)C(=O)O2)C(=O)N1")[0]) # True
print(is_penicillin("CC1(C(N2C(S1)C(C2=O)NC(=O)[C@H](N=C)C3=CC=CC=C3)C(=O)O)C")[0]) # True
print(is_penicillin("CC1(C(N2C(S1)C(C2=O)NC(=O)C(CC3=CC=CC=C3)OC)C(=O)O)C")[0]) # True