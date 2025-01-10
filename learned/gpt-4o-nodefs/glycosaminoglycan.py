"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is classified as a glycosaminoglycan based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool, str: True if molecule is likely a glycosaminoglycan, False otherwise
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for significant number of nitrogen atoms (amino component)
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if nitrogen_count < 2:
        return False, "Insufficient nitrogen atoms for typical glycosaminoglycan structure."

    # Check for a significant number of oxygen atoms (suggesting polysaccharides)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 8:
        return False, "Insufficient oxygen atoms for typical polysaccharide structure."

    # Check for sulfur atoms (sulfate groups)
    sulfur_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)

    # Since GAGs generally are long chains, check the number of rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Insufficient flexibility for a polymer structure"

    # Molecular weight typically higher due to the polysaccharide chain
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for a glycosaminoglycan."

    reason = "The molecule has characteristics consistent with glycosaminoglycans: "
    reason += f"Nitrogen atoms: {nitrogen_count}, Oxygen atoms: {oxygen_count}, "
    reason += f"Sulfur atoms: {sulfur_count}, Molecular weight: {mol_wt:.2f}"

    return True, reason

# Note: This is still a heuristic approach and might not cover all specific characteristics 
# and structural diversity of glycosaminoglycans.