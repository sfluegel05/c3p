"""
Classifies: CHEBI:138675 gas molecular entity
"""
"""
Classifies: CHEBI:138675 gas molecular entity
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a gas at STP based on SMILES string.
    Criteria: main group elements only, molecular weight <= 300, and absence of certain functional groups.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"
    
    # Define allowed main group elements (H, He, groups 1,2,13-18 up to Rn)
    allowed_atomic_numbers = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                              11, 12, 13, 14, 15, 16, 17, 18,
                              19, 20, 31, 32, 33, 34, 35, 36,
                              37, 38, 49, 50, 51, 52, 53, 54,
                              55, 56, 81, 82, 83, 84, 85, 86}
    
    # Check all atoms are main group
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_numbers:
            return False, f"Contains non-main group element {atom.GetSymbol()}"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 300:
        return False, f"Molecular weight {mol_wt:.2f} exceeds 300"
    
    # Check for functional groups that typically indicate non-gaseous state
    hydroxyl = Chem.MolFromSmarts("[OH]")
    carboxylic_acid = Chem.MolFromSmarts("C(=O)O")
    sulfonic_acid = Chem.MolFromSmarts("S(=O)(=O)O")
    
    if mol.HasSubstructMatch(hydroxyl):
        return False, "Contains hydroxyl group"
    if mol.HasSubstructMatch(carboxylic_acid):
        return False, "Contains carboxylic acid group"
    if mol.HasSubstructMatch(sulfonic_acid):
        return False, "Contains sulfonic acid group"
    
    # Additional check for ionic compounds (charges)
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0:
            return False, "Molecule contains charged species"
    
    return True, "Main group elements, MW <=300, no problematic groups"