"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: CHEBI:37405 phenyl acetates
An acetate ester obtained by formal condensation of the carboxy group of acetic acid
with the hydroxy group of any phenol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenyl acetate, False otherwise
        str: Reason for classification
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Look for aromatic ring(s)
        aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() if mol.GetRingInfo().IsAromaticRing(ring)]
        if not aromatic_rings:
            return False, "No aromatic rings found"

        # Look for acetate group(s) derived from acetic acid
        acetic_acid_pattern = Chem.MolFromSmarts("CC(=O)O")
        acetic_acid_matches = mol.GetSubstructMatches(acetic_acid_pattern)
        if not acetic_acid_matches:
            return False, "No acetate groups derived from acetic acid"

        # Check if any acetate group is attached to an aromatic ring (directly or through a linker)
        for acetate_match in acetic_acid_matches:
            acetate_oxygen = mol.GetAtomWithIdx(acetate_match[1])
            for neighbor in acetate_oxygen.GetNeighbors():
                if any(neighbor.IsInRing() for ring in aromatic_rings):
                    return True, "Contains a phenol ring with an acetate group derived from acetic acid attached"

        return False, "Acetate group not attached to an aromatic ring"

    except Exception as e:
        return False, f"Error occurred: {str(e)}"