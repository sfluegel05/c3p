"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: CHEBI:38170 quinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    A quinone is defined as a compound having a fully conjugated cyclic dione structure,
    such as that of benzoquinones, derived from aromatic compounds by conversion of an 
    even number of -CH= groups into -C(=O)- groups with any necessary rearrangement of 
    double bonds (polycyclic and heterocyclic analogues are included).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define quinone pattern
    quinone_pattern = Chem.MolFromSmarts("[#6]1([#6]=[#6])([#8])=[#6][#6]=[#6][#6]1=[#8]")
    
    # Check for quinone pattern
    if not mol.HasSubstructMatch(quinone_pattern):
        return False, "No quinone substructure found"
    
    # Check for conjugated cyclic dione
    diones = [m for m in mol.GetSubstructMatches(quinone_pattern)]
    is_conjugated = all(check_conjugation(mol, dione) for dione in diones)
    
    if not is_conjugated:
        return False, "Quinone substructure is not fully conjugated"
    
    # Check for ring membership
    ring_info = mol.GetRingInfo()
    is_cyclic = all(ring_info.IsAtomInRingOfSize(atom_idx, 6) for atom_idx in diones[0])
    
    if not is_cyclic:
        return False, "Quinone substructure is not cyclic"
    
    return True, "Molecule contains a fully conjugated cyclic dione structure (quinone)"

def check_conjugation(mol, dione):
    """
    Checks if a given dione substructure is fully conjugated in a molecule.
    
    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object
        dione (list): List of atom indices representing the dione substructure
        
    Returns:
        bool: True if the dione substructure is fully conjugated, False otherwise
    """
    for i in range(len(dione)):
        atom1 = mol.GetAtomWithIdx(dione[i])
        atom2 = mol.GetAtomWithIdx(dione[(i+1) % len(dione)])
        bond = mol.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx())
        if bond is None or bond.GetBondType() != Chem.BondType.DOUBLE:
            return False
    return True