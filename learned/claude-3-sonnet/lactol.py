"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies chemical entities of the class lactol:
Cyclic hemiacetals formed by intramolecular addition of a hydroxy group to an aldehydic or ketonic carbonyl group.
They are thus 1-oxacycloalkan-2-ols or unsaturated analogues.
"""

from rdkit import Chem
from rdkit.Chem import rdchem, rdMolDescriptors

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cyclic ethers
    cyclic_ethers = rdchem.GetSSSR(mol)
    if not cyclic_ethers:
        return False, "No cyclic ether found"

    # Check for hydroxyl groups adjacent to cyclic ethers
    for ring in cyclic_ethers:
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 8 and atom.GetDegree() == 2:  # Cyclic ether oxygen
                neighbors = [mol.GetAtomWithIdx(nbr_idx) for nbr_idx in atom.GetNeighbors()]
                if any(nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() == 1 for nbr in neighbors):
                    # Adjacent hydroxyl group found
                    return True, "Molecule contains a lactol structure"

    return False, "No lactol structure found"