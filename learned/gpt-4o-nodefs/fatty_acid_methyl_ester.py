"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for methyl ester group pattern (O=C(OC))
    ester_pattern = Chem.MolFromSmarts("O=C(OC)")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No methyl ester group found"
    
    # Count carbon atoms, excluding the methyl group in the ester
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6) - 1
    
    # Minimum number of carbon atoms typically found in fatty acids is 6 (excluding methyl for typical structures)
    if carbon_count < 6:
        return False, f"Chain too short for typical fatty acid methyl ester: {carbon_count + 1} carbons"
    
    # Check for presence of unwanted functional groups or complexity
    # Avoid structures with too many rings, multiple heteroatoms, or atypical saturation patterns
    # We can use RDKit's functions to filter out unlikely structures
    # Check for heavy atoms, there should ideally be just C, H, O
    if any(atom.GetAtomicNum() not in {1, 6, 8} for atom in mol.GetAtoms()):
        return False, "Contains atypical atoms for FAME"

    # Check number of rings as FAMEs are typically not cyclic
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return False, "Contains rings, which is unusual for a typical FAME"

    # FAMEs typically don't have exotic components
    return True, "Contains a methyl ester group and a sufficiently long typical carbon chain representing FAME characteristics"

__metadata__ = {
    'chemical_class': {
        'id': 'None',
        'name': 'fatty acid methyl ester',
        'definition': 'Esters derived from fatty acids and methanol, characterized by a methoxy group attached to an ester carbonyl.'
    }
}