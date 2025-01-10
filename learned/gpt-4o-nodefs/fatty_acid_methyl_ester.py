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
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No methyl ester group found"
    
    # Count carbon atoms, excluding the methyl group in the ester
    # Here, specifically target primary ester at the end, assuming methoxy -> carbon root
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6) - 1
    
    # Minimum number of carbon atoms typically found in fatty acids is 4 (lowered threshold for diversity)
    if carbon_count < 4:
        return False, f"Chain too short for typical fatty acid methyl ester: {carbon_count + 1} carbons"
    
    # Adjusted checks for complex structures
    # Check atoms and permit them with consideration to typical FAME variants
    atypical_atoms_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in {1, 6, 8})
    if atypical_atoms_count > 2:  # Allowing some variance
        return False, "Contains an excessive number of atypical atoms for FAME"

    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 2:  # Permitting some rings based on provided data
        return False, "Too many rings, which is unusual for a typical FAME"

    return True, "Contains a methyl ester group and a carbon chain typical for FAME characteristics"

__metadata__ = {
    'chemical_class': {
        'id': 'None',
        'name': 'fatty acid methyl ester',
        'definition': 'Esters derived from fatty acids and methanol, characterized by a methoxy group attached to an ester carbonyl.'
    }
}