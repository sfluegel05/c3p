"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: Epoxy Fatty Acid
"""
from rdkit import Chem

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid is a heterocyclic fatty acid containing an epoxide ring as part of its structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    ca_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not ca_matches:
        return False, "No carboxylic acid group found"

    # Check for epoxide ring (three-membered ring with one oxygen and two carbons)
    # This method searches for oxygen atoms in 3-membered rings
    ring_info = mol.GetRingInfo()
    epoxide_found = False

    for ring in ring_info.AtomRings():
        if len(ring) == 3:
            num_oxygen = 0
            num_carbon = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                atomic_num = atom.GetAtomicNum()
                if atomic_num == 8:
                    num_oxygen += 1
                elif atomic_num == 6:
                    num_carbon += 1
            if num_oxygen == 1 and num_carbon == 2:
                epoxide_found = True
                break

    if not epoxide_found:
        return False, "No epoxide ring found"

    # Optional: Check if the epoxide ring is connected to the fatty acid chain
    # Since fatty acids are long hydrocarbon chains, we can check the proportion of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    num_heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 1)

    if num_carbons < 12:
        return False, f"Molecule has only {num_carbons} carbons, not sufficient for a fatty acid"

    return True, "Molecule contains a carboxylic acid group and an epoxide ring"