"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
"""
Classifies: CHEBI:36347 2-enoyl-CoA
An unsaturated fatty acyl-CoA in which the S-acyl group contains a double bond between positions 2 and 3.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove explicit hydrogens and ignore stereochemistry
    mol = Chem.RemoveHs(mol)
    
    # Split molecule into acyl group and CoA group
    acyl_group, coa_group = split_acyl_coa(mol)
    
    # Check for CoA group
    if coa_group is None:
        return False, "Missing CoA group"
    
    # Check for unsaturation in acyl group
    if not has_unsaturation(acyl_group):
        return False, "Acyl group is saturated"
    
    # Check for double bond between positions 2 and 3
    if not has_double_bond_2_3(acyl_group):
        return False, "Double bond not between positions 2 and 3"
    
    # Check for fatty acid chain
    if not has_fatty_acid_chain(acyl_group):
        return False, "Missing fatty acid chain"
    
    # Check molecular weight (typically 800-1000 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800 or mol_wt > 1000:
        return False, "Molecular weight outside typical range"
    
    # Check atom counts
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    
    if c_count < 25 or o_count < 10 or n_count < 3 or p_count < 2 or s_count < 1:
        return False, "Atom counts outside typical range"
    
    return True, "Molecule is a 2-enoyl-CoA"

def split_acyl_coa(mol):
    """
    Splits a molecule into the acyl group and CoA group.

    Args:
        mol (Mol): RDKit mol object

    Returns:
        Mol, Mol: Acyl group and CoA group (or None if not found)
    """
    # Define SMARTS pattern for CoA group
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(=O)([O-])OP(=O)([O-])OC1C(OP(=O)([O-])[O-])C(OC2N=CN=C2N)C(OP(=O)([O-])[O-])C1N)C(=O)NCCC(=O)NCCS")
    
    # Match pattern and split molecule
    match = mol.GetSubstructMatches(coa_pattern)
    if not match:
        return mol, None
    
    acyl_atoms = list(set(range(mol.GetNumAtoms())) - set(match[0]))
    acyl_group = Chem.EditableMol(mol).GetMol()
    acyl_group = Chem.EditableMol(acyl_group).RemoveAtoms([atom for atom in range(mol.GetNumAtoms()) if atom not in acyl_atoms])
    acyl_group = acyl_group.GetMol()
    
    coa_atoms = match[0]
    coa_group = Chem.EditableMol(mol).GetMol()
    coa_group = Chem.EditableMol(coa_group).RemoveAtoms([atom for atom in range(mol.GetNumAtoms()) if atom not in coa_atoms])
    coa_group = coa_group.GetMol()
    
    return acyl_group, coa_group

def has_unsaturation(mol):
    """
    Checks if a molecule has at least one double bond.

    Args:
        mol (Mol): RDKit mol object

    Returns:
        bool: True if molecule has at least one double bond, False otherwise
    """
    adj_matrix = Chem.GetAdjacencyMatrix(mol)
    bonds = [sum(row) for row in adj_matrix]
    return max(bonds) > 3

def has_double_bond_2_3(mol):
    """
    Checks if a molecule has a double bond between positions 2 and 3 of the acyl chain.

    Args:
        mol (Mol): RDKit mol object

    Returns:
        bool: True if molecule has a double bond between positions 2 and 3, False otherwise
    """
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            
            # Check if atoms are part of the acyl chain
            if atom1.GetDegree() > 1 and atom2.GetDegree() > 1:
                # Get neighbor atoms
                neighbors1 = [atom.GetIdx() for atom in atom1.GetNeighbors()]
                neighbors2 = [atom.GetIdx() for atom in atom2.GetNeighbors()]
                
                # Check if neighbors are separated by 2 bonds
                if any(abs(n1 - n2) == 2 for n1 in neighbors1 for n2 in neighbors2):
                    return True
    
    return False

def has_fatty_acid_chain(mol):
    """
    Checks if a molecule has a fatty acid chain attached to the acyl group.

    Args:
        mol (Mol): RDKit mol object

    Returns:
        bool: True if molecule has a fatty acid chain, False otherwise
    """
    # Define SMARTS pattern for fatty acid chain
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    
    # Match pattern
    match = mol.GetSubstructMatches(fatty_acid_pattern)
    
    return len(match) > 0