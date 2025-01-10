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
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for methyl ester group pattern (O=C(OC))
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No methyl ester group found"

    # Check for a long hydrocarbon chain
    num_carbon_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbon_atoms < 8:  # Arbitrary choice for "long chain", can be adjusted
        return False, f"Chain too short for typical fatty acid methyl ester: {num_carbon_atoms} carbons"

    return True, "Contains a methyl ester group with a sufficiently long carbon chain"

__metadata__ = {
    'chemical_class': {
        'id': 'None',
        'name': 'fatty acid methyl ester',
        'definition': 'Esters derived from fatty acids and methanol, characterized by a methoxy group attached to an ester carbonyl.'
    }
}