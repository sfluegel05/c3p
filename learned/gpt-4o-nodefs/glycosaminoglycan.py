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
        bool, str: True if the molecule is likely a glycosaminoglycan, False otherwise, with a reason
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count nitrogen atoms (typical for amino components) and oxygen atoms (typical for polysaccharides)
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    sulfur_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)

    # Molecular weight calculation
    mol_wt = Descriptors.ExactMolWt(mol)

    # Check for number of rotatable bonds indicative of a polymer structure
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)

    # Define thresholds - these may need to be adjusted based on domain knowledge and datasets
    if nitrogen_count < 2:
        return False, "Insufficient nitrogen atoms for typical glycosaminoglycan structure."
    if oxygen_count < 10:
        return False, "Insufficient oxygen atoms for typical polysaccharide structure."
    if sulfur_count < 1:  # Assuming at least one sulfate group might be expected
        return False, "Insufficient sulfur atoms for typical sulfate groups in glycosaminoglycans."
    if n_rotatable < 15:
        return False, "Insufficient flexibility suggesting it's not a long polymer structure."
    if mol_wt < 500:  # Expecting higher molecular weight due to polysaccharides
        return False, "Molecular weight too low for a glycosaminoglycan."

    reason = "The molecule has characteristics consistent with glycosaminoglycans: "
    reason += f"Nitrogen atoms: {nitrogen_count}, Oxygen atoms: {oxygen_count}, "
    reason += f"Sulfur atoms: {sulfur_count}, Molecular weight: {mol_wt:.2f}, "
    reason += f"Rotatable bonds: {n_rotatable}"

    return True, reason

# Note: This is a heuristic approach and might need further refinement.