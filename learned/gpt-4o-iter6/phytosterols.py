"""
Classifies: CHEBI:26125 phytosterols
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols similar to cholesterol which occur in plants 
    and may vary in carbon side chains and/or presence or absence of a double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Steroid backbone pattern: More flexible tetracyclic core
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C1CCC3C2CCC4C3CCCC4")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No flexible tetracyclic steroid backbone found"
    
    # Phytosterols often have alcohol groups; checking for hydroxyl groups on the ring system
    if not any(atom.GetAtomicNum() == 8 and atom.GetTotalDegree() > 1 for atom in mol.GetAtoms()):
        return False, "No hydroxyl group found, but phytosterols usually have at least one"
    
    # Allow for presence or absence of double bonds anywhere in the structure
    # Phytosterols frequently have at least some unsaturation in the side chain or ring system
    unsaturation_patterns = [
        Chem.MolFromSmarts("C=C"),  # Simple olefin
        Chem.MolFromSmarts("C#C"),  # Acetylenic link not particularly common but included for extensiveness
    ]
    has_unsaturation = any(mol.HasSubstructMatch(pattern) for pattern in unsaturation_patterns)
    
    # Conclude based on presence of unsaturation and specific backbone checks
    if has_unsaturation:
        return True, "Contains elements of unsaturation and flexible tetracyclic steroid backbone typical of phytosterols"
    
    return False, "Fails to meet unsaturation and backbone criteria typical of phytosterols"