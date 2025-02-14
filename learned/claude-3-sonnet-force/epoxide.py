"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: CHEBI:33856 epoxide

An epoxide is any cyclic ether in which the oxygen atom forms part of a 3-membered ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxide, False otherwise
        str: Reason for classification
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Look for epoxide ring patterns
        epoxide_patterns = [
            Chem.MolFromSmarts("[O;R]1[CR]2[CR]1[CR]2"),  # Basic epoxide ring
            Chem.MolFromSmarts("[O;R]1[CR]2[NR]1[CR]2"),  # Epoxide ring with N in the ring
            Chem.MolFromSmarts("[O;R]1[CR]2[OR]1[CR]2"),  # Epoxide ring with O in the ring
            Chem.MolFromSmarts("[O;R]1[CR]2[SR]1[CR]2"),  # Epoxide ring with S in the ring
        ]

        # Check for epoxide ring matches
        has_epoxide = any(mol.HasSubstructMatch(pattern) for pattern in epoxide_patterns)
        if not has_epoxide:
            return False, "No epoxide ring found"

        # Additional checks to filter out false positives
        epoxide_atoms = [atom for atom in mol.GetAtoms() if atom.IsInRingSize(3) and atom.GetAtomicNum() == 8]
        if not epoxide_atoms:
            return False, "No epoxide oxygen atom found"

        for epoxide_atom in epoxide_atoms:
            neighbors = [mol.GetAtomWithIdx(n).GetAtomicNum() for n in epoxide_atom.GetNeighbors()]
            if 6 not in neighbors or len(neighbors) != 2:
                return False, "Epoxide oxygen not connected to two carbon atoms"

        # Count the number of epoxide rings
        epoxide_rings = mol.GetRingInfo().AtomRings()
        epoxide_rings = [ring for ring in epoxide_rings if len(ring) == 3 and mol.GetAtomWithIdx(ring[0]).GetAtomicNum() == 8]
        n_epoxide_rings = len(epoxide_rings)

        if n_epoxide_rings > 0:
            return True, f"Contains {n_epoxide_rings} epoxide ring(s)"
        else:
            return False, "No valid epoxide ring found"

    except Exception as e:
        return False, f"Error: {str(e)}"