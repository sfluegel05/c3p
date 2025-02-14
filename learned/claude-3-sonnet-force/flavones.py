"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: CHEBI:17794 flavones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    A flavone is a flavonoid with a 2-aryl-1-benzopyran-4-one (2-arylchromen-4-one) skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for flavone skeleton pattern and its tautomers
    flavone_patterns = [
        Chem.MolFromSmarts("c1cc(oc2ccccc2c1=O)-c"),  # 2-aryl-1-benzopyran-4-one
        Chem.MolFromSmarts("c1cc(oc2ccccc2C(=O)O)-c")  # Tautomer with OH group
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in flavone_patterns):
        return False, "Does not contain the flavone skeleton or its tautomers"
    
    # Count aromatic rings
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if aromatic_rings != 3:
        return False, f"Expected 3 aromatic rings, found {aromatic_rings}"
    
    # Check for oxygens
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygens < 4:
        return False, "Too few oxygens for a flavone"
    
    # Check for hydroxyl groups
    hydroxy_pattern = Chem.MolFromSmarts("OC")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    num_hydroxy = len(hydroxy_matches)
    if num_hydroxy < 1:
        return False, "No hydroxyl groups found"
    
    # Check for common flavone substituents and glycosidic substituents
    substituents = ["CH3", "OCH3", "OH", "O", "CH2", "CH", "OC1C(O)C(O)C(O)C(O)C1O"]  # Glucose
    glycosidic_pattern = Chem.MolFromSmarts("OC1OC(CO)C(O)C(O)C1O")  # Rhamnose
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if glycosidic_matches:
        substituents.append("OC1OC(CO)C(O)C(O)C1O")  # Add rhamnose if found
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic():
            neighbors = [nbr.GetAtomicNum() for nbr in atom.GetNeighbors()]
            substituent = "".join([get_symbol(n) for n in sorted(neighbors)])
            if substituent not in substituents:
                return False, f"Unusual substituent '{substituent}' found"
    
    # Check molecular weight and Lipinski's rules (optional)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 600:
        return False, "Molecular weight outside the typical range for flavones"
    
    if not rdMolDescriptors.CalcLipinski(mol):
        return False, "Violates Lipinski's rules for drug-likeness"
    
    return True, "Contains the flavone skeleton with expected aromatic rings, oxygens, substituents, and properties"

def get_symbol(atomic_num):
    return Chem.Atom(atomic_num).GetSymbol()