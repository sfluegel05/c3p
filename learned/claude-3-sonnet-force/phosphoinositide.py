"""
Classifies: CHEBI:18179 phosphoinositide
"""
"""
Classifies: CHEBI:17855 phosphoinositide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol that is phosphorylated at one or more of the hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for phosphatidylinositol backbone
    pi_pattern = Chem.MolFromSmarts("C(OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O)OC)=O")
    if not mol.HasSubstructMatch(pi_pattern):
        return False, "No phosphatidylinositol backbone found"
    
    # Check for two fatty acid chains
    fatty_acid_chains = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "O" and atom.GetDegree() == 2:
            neighbor1, neighbor2 = atom.GetNeighbors()
            if neighbor1.GetSymbol() == "C" and neighbor1.GetDegree() == 3 and neighbor2.GetSymbol() == "C" and neighbor2.GetDegree() == 2:
                chain = [neighbor1.GetIdx()]
                curr_atom = neighbor1
                while curr_atom.GetDegree() > 1:
                    neighbors = [n for n in curr_atom.GetNeighbors() if n.GetIdx() != chain[-1]]
                    if len(neighbors) == 0:
                        break
                    next_atom = neighbors[0]
                    chain.append(next_atom.GetIdx())
                    curr_atom = next_atom
                if len(chain) >= 6:  # Minimum chain length of 6 carbon atoms
                    fatty_acid_chains.append(chain)
    if len(fatty_acid_chains) != 2:
        return False, f"Found {len(fatty_acid_chains)} fatty acid chains, expected 2"
    
    # Check for unsaturation in fatty acid chains
    unsaturated_chains = []
    for chain in fatty_acid_chains:
        unsaturated = False
        for i in range(len(chain) - 1):
            atom1 = mol.GetAtomWithIdx(chain[i])
            atom2 = mol.GetAtomWithIdx(chain[i + 1])
            bond = mol.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx())
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                unsaturated = True
                break
        if unsaturated:
            unsaturated_chains.append(chain)
    if len(unsaturated_chains) == 0:
        return False, "No unsaturated fatty acid chains found"
    
    # Check for phosphate groups on the inositol ring
    inositol_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.IsInRingSize(6)]
    phosphate_groups = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == "P" and any(n.GetIdx() in inositol_atoms for n in atom.GetNeighbors())]
    if not phosphate_groups:
        return False, "No phosphate groups on the inositol ring"
    
    return True, f"Contains phosphatidylinositol backbone with {len(phosphate_groups)} phosphate group(s) on the inositol ring"