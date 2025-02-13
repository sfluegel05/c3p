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
    ester_pattern = Chem.MolFromSmarts("O=C(OC)")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No methyl ester group found"

    # Count carbon atoms, excluding the methyl group in the ester
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6) - 1

    # Minimum number of carbon atoms typically found in fatty acids is 4 (excluding methyl)
    if carbon_count < 4:
        return False, f"Chain too short for typical fatty acid methyl ester: {carbon_count + 1} carbons"

    # Check for presence of unwanted functional groups:
    # Add logic here to exclude molecules exhibiting complex structures atypical for simple FAMEs.
    # This might include looking for heteroatoms or unsaturations that aren't characteristic of typical fatty acid chains.

    return True, "Contains a methyl ester group with a sufficiently long carbon chain"

__metadata__ = {
    'chemical_class': {
        'id': 'None',
        'name': 'fatty acid methyl ester',
        'definition': 'Esters derived from fatty acids and methanol, characterized by a methoxy group attached to an ester carbonyl.'
    }
}