"""
Classifies: CHEBI:37581 gamma-lactone
"""
"""
Classifies: CHEBI:35524 gamma-lactone
A lactone having a five-membered lactone ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is a lactone with a five-membered ring containing an unsaturated bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gamma-lactone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for gamma-lactone pattern
    gamma_lactone_pattern = Chem.MolFromSmarts("[O;R]1[C;R]=C[C;R]=[C;R][C;R]=O1")
    gamma_lactone_matches = mol.GetSubstructMatches(gamma_lactone_pattern)
    
    if len(gamma_lactone_matches) == 0:
        return False, "No gamma-lactone substructure found"
    
    # Check if the matched substructure is a valid gamma-lactone
    for match in gamma_lactone_matches:
        ring_atoms = mol.GetAtomRingInfo().AtomRings()[match]
        if len(ring_atoms) != 5:
            continue
        
        # Check for ester group (-O-C=O-) in the ring
        ester_group_found = False
        for i in range(len(ring_atoms)):
            atom1 = mol.GetAtomWithIdx(ring_atoms[i])
            atom2 = mol.GetAtomWithIdx(ring_atoms[(i+1)%5])
            atom3 = mol.GetAtomWithIdx(ring_atoms[(i+2)%5])
            if (atom1.GetSymbol() == 'O' and atom2.GetSymbol() == 'C' and
                atom3.GetSymbol() == 'O' and atom2.GetFormalCharge() == 0):
                ester_group_found = True
                break
        
        if not ester_group_found:
            continue
        
        # Check for unsaturated bond between carbon atoms in the ring
        unsaturated_bond = False
        for i in range(len(ring_atoms)):
            atom1 = mol.GetAtomWithIdx(ring_atoms[i])
            atom2 = mol.GetAtomWithIdx(ring_atoms[(i+1)%5])
            if atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'C':
                bond = mol.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx())
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    unsaturated_bond = True
                    break
        
        if unsaturated_bond:
            return True, "Contains a five-membered lactone ring with an unsaturated bond"
    
    return False, "No valid gamma-lactone substructure found"