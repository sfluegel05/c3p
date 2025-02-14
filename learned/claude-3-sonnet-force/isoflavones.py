"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: CHEBI:28820 isoflavone

Any isoflavonoid with a 3-aryl-1-benzopyran-4-one (3-aryl-4H-chromen-4-one) skeleton and its substituted derivatives.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_isoflavone(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for isoflavone skeleton pattern (3-aryl-1-benzopyran-4-one)
    isoflavone_pattern = Chem.MolFromSmarts("c1cc(-c2oc3c(c2=O)c(O)cc(O)c3)ccc1")
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavone skeleton found"
    
    # Check substituents
    allowed_substituents = ['O', 'C', 'OC', 'Br', 'Cl', 'S', 'OS', 'CO', 'COC', 'OCC', 'OCCO', 'COCC', 'OCCOC', 'COCO']
    substituents = [atom.GetSmarts() for atom in mol.GetAtoms() if atom.GetSymbol() not in ['H', 'C', 'O'] and atom.GetSmarts() not in allowed_substituents]
    if substituents:
        return False, f"Found disallowed substituents: {', '.join(substituents)}"
    
    # Check molecular weight - isoflavones typically 250-300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 500:
        return False, "Molecular weight outside typical range for isoflavones"
    
    # Count oxygens, carbons, and heteroatoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    hetero_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6, 8])
    
    if o_count < 3 or o_count > 10:
        return False, "Oxygen count outside typical range for isoflavones"
    if c_count < 15 or c_count > 30:
        return False, "Carbon count outside typical range for isoflavones"
    if hetero_count > 5:
        return False, "Too many heteroatoms for isoflavone"
    
    return True, "Molecule contains the isoflavone skeleton with allowed substituents"