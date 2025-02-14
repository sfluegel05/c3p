"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
"""
Classifies: CHEBI:38698 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom, rdMolDescriptors

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    A 17alpha-hydroxy steroid is a steroid with a hydroxyl group attached at the 17th position
    in the alpha configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid scaffold
    steroid_pattern = Chem.MolFromSmarts("[C@]123[C@H]([C@@H]4[C@@]([C@@H]([C@H]([C@@H]2C)C(C)(C)C)C)(C)C=C5C[C@]6([H])C[C@@]7([H])C[C@@H](O)[C@]46C(=O)C=C7[C@]35C)O")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid scaffold found"
    
    # Embed molecule in 3D space
    AllChem.EmbedMolecule(mol)
    
    # Find the 17-hydroxyl group
    hydroxyl_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1]
    if len(hydroxyl_atoms) == 0:
        return False, "No hydroxyl group found"
    
    # Check alpha stereochemistry at C17
    for hydroxyl_idx in hydroxyl_atoms:
        atom = mol.GetAtomWithIdx(hydroxyl_idx)
        neighbors = [mol.GetAtomWithIdx(n).GetIdx() for n in atom.GetNeighbors()]
        if len(neighbors) == 1:
            c17_idx = neighbors[0]
            c17_atom = mol.GetAtomWithIdx(c17_idx)
            if c17_atom.GetHybridization() == Chem.HybridizationType.SP3:
                perm = list(range(mol.GetNumAtoms()))
                stereo = rdMolDescriptors.GetStereoWiggleCalculatorForConfId(mol).calculateStereoWiggle(c17_idx, perm)
                if stereo == 1:
                    break
    else:
        return False, "Hydroxyl group not in alpha configuration at C17"
    
    # Additional checks for steroid-like properties
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 4:
        return False, "Does not have 4 rings"
    
    if sum(1 for ring in ring_info.AtomRings() if ring.IsAromatic()) > 0:
        return False, "Contains aromatic rings, steroids should be fully saturated"
    
    return True, "Steroid with hydroxyl group in alpha configuration at C17"