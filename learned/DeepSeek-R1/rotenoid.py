"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: CHEBI:72564 rotenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    Rotenoids have a cis-fused tetrahydrochromeno[3,4-b]chromene skeleton with possible substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a rotenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Core scaffold pattern for cis-fused tetrahydrochromeno[3,4-b]chromene skeleton
    # Matches fused bicyclic system with oxygen bridges and ketone group
    scaffold_pattern = Chem.MolFromSmarts("[#6]1-&@[#8]-&@[#6]-&@[#6]2-&@[#6](=O)-&@[#6]-&@[#6]3-&@[#6]-&@[#6]-&@[#6]-&@[#6]-&@[#6]-3-&@[#6]-2-&@[#6]-1")
    if not mol.HasSubstructMatch(scaffold_pattern):
        return False, "Missing core chromenochromene scaffold"
    
    # Verify presence of at least one aromatic ring (characteristic of chromene systems)
    aromatic_rings = [ring for ring in Chem.GetSymmSSSR(mol) if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
    if len(aromatic_rings) < 2:
        return False, "Insufficient aromatic rings for chromene systems"
    
    # Check for oxygen atoms in rings (characteristic of chromene oxygen bridges)
    oxygen_in_rings = any(any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring) for ring in Chem.GetSymmSSSR(mol))
    if not oxygen_in_rings:
        return False, "Missing oxygen in ring systems"
    
    return True, "Contains cis-fused chromenochromene scaffold with oxygen bridges"