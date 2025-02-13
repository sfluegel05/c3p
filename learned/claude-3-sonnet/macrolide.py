"""
Classifies: CHEBI:25106 macrolide
"""
"""
Classifies: CHEBI:37219 macrolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide is a macrocyclic lactone with a ring of twelve or more members derived from a polyketide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for macrocyclic ring with >=12 atoms
    ring_info = mol.GetRingInfo()
    has_macrocycle = any(len(ring) >= 12 for ring in ring_info.AtomRings())
    if not has_macrocycle:
        return False, "No macrocyclic ring with 12 or more atoms found"
    
    # Check for lactone group(s)
    lactone_patterns = ['[C&r]1[C&r]([C&r](=[O&r])[O&r]1)', # Simple lactone
                        '[C&r]1[C&r]([C&r](=[O&r])[O&r][C&r]1)', # Lactone with ring substituent
                        '[C&r]1[C&r]2[C&r]([C&r]([C&r]2=[O&r])[O&r]1)']  # Spiro lactone
    has_lactone = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in lactone_patterns)
    if not has_lactone:
        return False, "No lactone group found"
    
    # Check for polyketide-derived structure
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_carbonyls = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and sum(mol.GetAtomWithIdx(i).GetTotalNumHs() for i in atom.GetNeighbors()) == 0)
    
    if n_carbons < 12 or n_oxygens < 3 or n_carbonyls < 2:
        return False, "Structure does not appear to be polyketide-derived"
    
    return True, "Contains a macrocyclic lactone ring derived from a polyketide"