"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol phosphate in which the inositol component has myo-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all rings of size 6 (potential cyclohexane rings)
    ssr = Chem.GetSymmSSSR(mol)
    cyclohexane_rings = [ring for ring in ssr if len(ring) == 6]
    if not cyclohexane_rings:
        return False, "No cyclohexane rings found"

    # For each cyclohexane ring, check if all carbons are connected to oxygen (possibly substituted)
    found_inositol_phosphate = False
    for ring in cyclohexane_rings:
        carbons_with_oxygen = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                break  # Not a carbon atom
            has_oxygen = False
            oxygen_atom = None
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    has_oxygen = True
                    oxygen_atom = neighbor
                    break
            if not has_oxygen:
                break  # This carbon doesn't have an attached oxygen
            carbons_with_oxygen += 1
        if carbons_with_oxygen == 6:
            found_inositol_phosphate = True
            break
    if not found_inositol_phosphate:
        return False, "No inositol-like ring found"

    # Optionally, check for phosphate groups attached to the oxygens
    has_phosphate = False
    for idx in ring:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                oxygen_atom = neighbor
                for neighbor2 in oxygen_atom.GetNeighbors():
                    if neighbor2.GetAtomicNum() == 15:  # Phosphorus
                        has_phosphate = True
                        break
                if has_phosphate:
                    break
        if has_phosphate:
            break
    if not has_phosphate:
        return False, "No phosphate groups attached to inositol ring"

    return True, "Contains inositol phosphate core structure"