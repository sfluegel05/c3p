"""
Classifies: CHEBI:74716 withanolide
"""
"""
Classifies: withanolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Fragments

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide is defined as any steroid lactone that is a C28 steroid 
    with a modified side chain forming a lactone ring and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if molecule has a steroid nucleus using RDKit fragment function
    if Fragments.fr_steroid(mol) == 0:
        return False, "No steroid nucleus found"
    
    # Define general lactone pattern (any cyclic ester)
    lactone_pattern = Chem.MolFromSmarts("*1CCOC(=O)C1")  # 5-membered lactone
    lactone_pattern6 = Chem.MolFromSmarts("*1CCCCOC(=O)C1")  # 7-membered lactone
    ester_in_ring = False
    ri = mol.GetRingInfo()
    for ring_atoms in ri.AtomRings():
        ring = Chem.PathToSubmol(mol, ring_atoms)
        # Check if ring contains an ester group
        if ring.HasSubstructMatch(Chem.MolFromSmarts("C(=O)O[C;R]")):
            ester_in_ring = True
            break
    if not ester_in_ring:
        return False, "No lactone ring found"
    
    # Check connection between steroid nucleus and lactone ring
    # Get the steroid nucleus atoms
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4C3(CCCC4)C')  # Steroid nucleus
    steroid_match = mol.GetSubstructMatch(steroid_pattern)
    if not steroid_match:
        return False, "No steroid nucleus found"
    steroid_atoms = set(steroid_match)
    
    # Find ester bonds in rings
    ester_bond = Chem.MolFromSmarts('C(=O)O[C;R]')
    ester_bonds = mol.GetSubstructMatches(ester_bond)
    lactone_in_side_chain = False
    for bond_match in ester_bonds:
        bond_atoms = set(bond_match)
        # Check if ester bond is in a ring
        atom1 = mol.GetAtomWithIdx(bond_match[0])
        atom2 = mol.GetAtomWithIdx(bond_match[1])
        if ri.IsAtomInRingOfSize(atom1.GetIdx(), 5) or ri.IsAtomInRingOfSize(atom1.GetIdx(), 6) or \
           ri.IsAtomInRingOfSize(atom1.GetIdx(), 7):
            # Check if ester is connected to steroid nucleus via side chain
            if not bond_atoms.issubset(steroid_atoms):
                lactone_in_side_chain = True
                break
    if not lactone_in_side_chain:
        return False, "Lactone ring is not in the side chain"
    
    # Optional: check carbon count, allow some flexibility
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 26:
        return False, f"Too few carbons for a withanolide (found {c_count} carbons)"
    
    return True, "Molecule is a withanolide (steroid nucleus with side-chain lactone ring)"