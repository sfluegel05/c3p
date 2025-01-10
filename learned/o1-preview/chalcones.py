"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: chalcones
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is a ketone that is 1,3-diphenylpropenone (benzylideneacetophenone), ArCH=CH(=O)Ar, and its derivatives formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define chalcone core pattern: aromatic ring connected via alpha,beta-unsaturated ketone to another aromatic ring
    chalcone_pattern = Chem.MolFromSmarts("[a][C]=[C]-C(=O)-[a]")
    if chalcone_pattern is None:
        return False, "Invalid chalcone SMARTS pattern"

    # Check for the chalcone core
    if not mol.HasSubstructMatch(chalcone_pattern):
        return False, "Chalcone core not found"

    # Optionally, verify the number of aromatic rings
    ri = mol.GetRingInfo()
    num_aromatic_rings = sum(1 for ring in ri.AtomRings() if all(mol.GetAtomWithIdx(atom_idx).GetIsAromatic() for atom_idx in ring))
    if num_aromatic_rings < 2:
        return False, "Less than two aromatic rings detected"

    return True, "Chalcone core detected"