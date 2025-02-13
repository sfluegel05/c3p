"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies chemical entities of the class lactol:
Cyclic hemiacetals formed by intramolecular addition of a hydroxy group to an aldehydic or ketonic carbonyl group.
They are thus 1-oxacycloalkan-2-ols or unsaturated analogues.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

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

    # Look for lactol pattern: cyclic ether with OH group on adjacent carbon
    lactol_pattern = Chem.MolFromSmarts("[OH,OH2][C;r]1[O;r][C;r]=[C;r]~[C;r]~[C;r]~[C;r]1")
    matches = mol.GetSubstructMatches(lactol_pattern)

    if matches:
        # Check if the matched substructure is a valid lactol
        for match in matches:
            atoms = [mol.GetAtomWithIdx(idx) for idx in match]
            ring_info = mol.GetRingInfo().AtomRings()

            # Check if the atoms form a ring and if the ring contains a cyclic ether
            if any(set(match) <= ring for ring in ring_info):
                has_ether = any(atom.GetAtomicNum() == 8 and atom.GetDegree() == 2 for atom in atoms)
                if has_ether:
                    return True, "Molecule contains a lactol structure"

    return False, "No lactol structure found"