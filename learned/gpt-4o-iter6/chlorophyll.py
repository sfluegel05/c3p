"""
Classifies: CHEBI:28966 chlorophyll
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    Chlorophylls are defined by a magnesium porphyrin core with a fifth ring, and usually a long phytol chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chlorophyll, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for magnesium atom
    magnesium_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'Mg']
    if len(magnesium_atoms) != 1:
        return False, "Magnesium not found or more than one magnesium atom present"

    # Check for porphyrin core (four pyrrole-like N-atoms)
    porphyrin_n_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'N']
    if len(porphyrin_n_atoms) < 4:
        return False, "Less than four nitrogen atoms found"

    # Check for the fifth ring
    # Look for a ring system size larger than 4
    ssr = Chem.GetSymmSSSR(mol)
    five_rings = 0
    for ring in ssr:
        if len(ring) > 4:
            five_rings += 1
    
    if five_rings == 0:
        return False, "Fifth ring not detected"

    # Check for the presence of a phytol chain (at least a long aliphatic chain)
    # Using a simple rotatable bond or carbon count could help detect this
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C']
    if len(carbon_atoms) < 20:  # adjust number depending on expected chain length
        return False, "Too few carbon atoms, long hydrophobic chain possibly absent"

    # Confirm long chain through rotatable bonds (for ample flexibility)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:  # adjust number depending on expected chain length
        return False, "Not enough rotatable bonds for a phytol chain"

    return True, "Contains magnesium chelated by four pyrrole-like rings, a fifth ring, and a long hydrophobic chain"