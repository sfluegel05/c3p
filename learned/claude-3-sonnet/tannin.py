"""
Classifies: CHEBI:26848 tannin
"""
"""
Classifies: CHEBI:18035 tannin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are polyphenolic compounds that typically contain galloyl or hexahydroxydiphenoyl moieties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for galloyl substructure
    galloyl_pattern = Chem.MolFromSmarts("O=C(O)c1cc(O)c(O)c(O)c1")
    if mol.HasSubstructMatch(galloyl_pattern):
        return True, "Contains galloyl moiety, characteristic of tannins"
    
    # Look for hexahydroxydiphenoyl substructure
    hhdp_pattern = Chem.MolFromSmarts("O=C(O)c1c(O)c(O)c(O)c(C(=O)O)c1O")
    if mol.HasSubstructMatch(hhdp_pattern):
        return True, "Contains hexahydroxydiphenoyl moiety, characteristic of tannins"
    
    # Check for high degree of hydroxylation and polyphenolic nature
    n_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    n_hydroxyl = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
    
    if n_aromatic_rings > 1 and n_hydroxyl >= 5:
        return True, "Polyphenolic compound with high degree of hydroxylation, likely a tannin"
    
    return False, "No characteristic tannin substructures found"