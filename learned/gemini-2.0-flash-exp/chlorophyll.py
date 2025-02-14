"""
Classifies: CHEBI:28966 chlorophyll
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    Chlorophylls are magnesium porphyrins with a fifth ring and usually a long phytol chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chlorophyll, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for magnesium
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 12]
    if not mg_atoms:
        return False, "No magnesium atom found."
    if len(mg_atoms) != 1:
        return False, "Incorrect number of magnesium atoms."

    # 2. Flexible porphyrin core check (4 nitrogens coordinating with Mg)
    porphyrin_core_pattern = Chem.MolFromSmarts("[nX2][c][c][n][c][c][n][c][c][n][Mg]")
    if not mol.HasSubstructMatch(porphyrin_core_pattern):
       return False, "No porphyrin core with magnesium found."

    # 3. Check for a fifth ring attached to the porphyrin core. We look for a fused ring connected to a carbon
    #    of the porphyrin core. This allows for 5- or 6-membered rings, or even other fused ring systems
    fifth_ring_pattern = Chem.MolFromSmarts("[c]12[c][c][c]([c][c]1)[c]2")
    if not mol.HasSubstructMatch(fifth_ring_pattern):
        return False, "No fifth ring attached to the porphyrin core found."

    # 4. Check for phytol or similar long chain (at least 1 ester and 1 long chain)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not ester_matches:
        return False, "No ester group found."
    
    # Phytol chain: long carbon chain with branches. We also accept generic long fatty chain
    phytol_pattern = Chem.MolFromSmarts("[CX4]([CX4])~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    long_chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")

    has_long_chain = False
    for match in ester_matches:
        ester_oxygen_atom = mol.GetAtomWithIdx(match[0])

        for neighbor in ester_oxygen_atom.GetNeighbors():
           if mol.HasSubstructMatch(phytol_pattern, fromAtomIdx = neighbor.GetIdx()) or mol.HasSubstructMatch(long_chain_pattern, fromAtomIdx = neighbor.GetIdx()):
             has_long_chain = True
             break
        if has_long_chain:
           break
           
    if not has_long_chain:
        return False, "Missing phytol or long chain."

    
    
    return True, "Contains a magnesium porphyrin core, a fifth ring, and a long chain"