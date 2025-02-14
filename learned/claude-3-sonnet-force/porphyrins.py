"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: CHEBI:28799 porphyrins
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdMolTransforms

def is_porphyrin(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    Porphyrins are defined as natural pigments containing a fundamental skeleton
    of four pyrrole nuclei united through the alpha-positions by four methine
    groups to form a macrocyclic structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for macrocyclic structure
    if not mol.getNumConformers():
        mol = Chem.AddHs(mol)
        mol = Chem.RemoveHs(mol)  # Add/remove hydrogens to force ring perception
    if not mol.getNumConformers():
        return False, "Not a cyclic structure"
    
    # Look for 4 pyrrole rings
    pyrrole_pattern = Chem.MolFromSmarts("[nH]1cccc1")
    pyrrole_rings = mol.GetSubstructMatches(pyrrole_pattern)
    if len(pyrrole_rings) != 4:
        return False, f"Found {len(pyrrole_rings)} pyrrole rings, need exactly 4"
    
    # Check for connection via methine bridges
    methine_bridge_pattern = Chem.MolFromSmarts("[CH]=")
    methine_bridges = mol.GetSubstructMatches(methine_bridge_pattern)
    if len(methine_bridges) != 4:
        return False, f"Found {len(methine_bridges)} methine bridges, need exactly 4"
    
    # Check planarity
    planar = rdMolDescriptors.CalcPlanarityFast(mol)
    if planar < 0.9:
        return False, "Structure is not planar"
    
    # Check for optional metal coordination
    metal_atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 20]
    if metal_atoms:
        if len(metal_atoms) > 1:
            return False, "Multiple metal atoms detected"
        metal_symbol = Chem.GetPeriodicTable().GetElementSymbol(metal_atoms[0])
        return True, f"Porphyrin macrocycle with coordinated {metal_symbol} atom"
    
    return True, "Porphyrin macrocycle"

# Example usage
print(is_porphyrin("Fc1c(F)c(F)c(c(F)c1F)C1=C2C=CC3=[N]2[Pd]24n5c1ccc5C(c1c(F)c(F)c(F)c(F)c1F)=C1C=CC(C(c5c(F)c(F)c(F)c(F)c5F)=c5ccc(=C3c3c(F)c(F)c(F)c(F)c3F)n25)=[N]41"))
# Output: (True, 'Porphyrin macrocycle with coordinated Pd atom')

print(is_porphyrin("C[C@]1(CC(O)=O)[C@H](CCC(O)=O)C2=CC3=[N+]4C(Cc5c(CCC(O)=O)c(CC(O)=O)c6C=C7[N+]8=C(C=C1N2[Co--]48n56)[C@@H](CCC(O)=O)[C@]7(C)CC(O)=O)=C(CCC(O)=O)C3CC(O)=O"))
# Output: (True, 'Porphyrin macrocycle with coordinated Co atom')