"""
Classifies: CHEBI:36092 clavulone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    Clavulones are a class of esterified prostanoids obtained from marine corals.
    They have a prostanoid core with several ester groups, a long chain, a hydroxyl
    group and halogens

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core prostanoid pattern (bicyclic ring with ketone on the 5-membered ring). Specific stereochemistry
    prostanoid_core_pattern = Chem.MolFromSmarts("[C]12[C@H]([C@@H]([C@@H]([C]1=O)C[C@H]2)([#6]))")
    if not mol.HasSubstructMatch(prostanoid_core_pattern):
        return False, "No prostanoid core found"

    # check that there are two rings
    ring_info = mol.GetRingInfo()
    if len(ring_info.AtomRings()) !=2:
        return False, "The molecule does not have the required two rings"


    # Check for any Oxygen atom in position 3, 8 or 12 of the prostanoid ring
    oxygen_on_prostanoid = False
    for match in mol.GetSubstructMatches(prostanoid_core_pattern):
       
        for atom_index in match:
            atom = mol.GetAtomWithIdx(atom_index)
            if atom.GetSymbol() == 'O':
                oxygen_on_prostanoid = True
                break
    
    if not oxygen_on_prostanoid:
        return False, "Prostanoid core must have an oxygen in any of the C-3, C-8, or C-12 positions"

    # Count ester groups (at least 2)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
          return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Check for long alkyl chain (at least C7)
    alkyl_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    alkyl_matches = mol.GetSubstructMatches(alkyl_chain_pattern)
    if len(alkyl_matches) == 0:
       return False, "No long alkyl chain detected"
   
    # Check for hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group present"
    
    # Check if there is any halogen
    halogen_present = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in [9, 17, 35, 53]:  # F, Cl, Br, I
             halogen_present = True
             break
    
    if not halogen_present:
       return False, "No halogen present, clavulones must have halogens"
    
    # Verify number of carbon atoms is greater than 15
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, "Too few carbon atoms for a clavulone"

    return True, "Matches clavulone characteristics (prostanoid core, ester groups, long chain, hydroxyl group and a halogen)"