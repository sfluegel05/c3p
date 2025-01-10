"""
Classifies: CHEBI:48039 dihydroflavonols
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol typically contains a characteristic flavanone core structure
    with hydroxylation and specific chiral centers.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the flavanone core structure pattern (3-ring system with chirality)
    flavanone_pattern = Chem.MolFromSmarts("C1([C@@H]2[C@H](C(=O)c3cc(O)ccc23)c2cc(O)ccc12)O")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone core with required chirality found"
    
    # Look for hydroxylation pattern on rings
    hydroxyl_benzene = Chem.MolFromSmarts("c1cc(O)c(O)c(O)c1")
    if not mol.HasSubstructMatch(hydroxyl_benzene):
        return False, "Hydroxylation pattern typical of dihydroflavonols not found"
    
    # Number of aromatic rings - at least 2 are needed
    n_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if n_aromatic_rings < 2:
        return False, "Fewer than 2 aromatic rings, typical of non-flavonoid structures"
    
    # Validate that it matches known hydroxylated dihydroflavonol structures
    n_hydroxyls = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1])
    if n_hydroxyls < 4:
        return False, "Insufficient hydroxyl groups, expected at least 4"
    
    return True, "Contains dihydroflavonol features including flavanone core with proper hydroxylation and chirality"