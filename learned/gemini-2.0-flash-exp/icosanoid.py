"""
Classifies: CHEBI:23899 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    Icosanoids are signalling molecules arising from the oxidation of C20 EFAs.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for C20 Carbon Backbone: Look for a chain of at least 18 C atoms to allow for branching and cyclic structures
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18:
        return False, f"Too few carbon atoms. Found: {c_count}, require at least 18"

    # 2. Check for at least 2 double bonds
    double_bond_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
    if double_bond_count < 2:
        return False, f"Too few double bonds. Found: {double_bond_count}, require at least 2"
    
    # 3. Check for Carboxylic Acid Group:
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    if not (mol.HasSubstructMatch(carboxylic_acid_pattern) or mol.HasSubstructMatch(ester_pattern)):
       return False, "No carboxylic acid or ester group found."

    # 4. Check for Oxygen Atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
      return False, f"Too few oxygen atoms. Found {o_count}, require at least 2."

    # 5 & 6. Specific Substructures and pattern matching.
    # Check for cyclopentane ring
    cyclopentane_pattern = Chem.MolFromSmarts("C1CCCC1")
    if mol.HasSubstructMatch(cyclopentane_pattern):
        
        # If it is a prostaglandin or thromboxane, it must have at least one -OH group on that ring.
        hydroxyl_on_ring = Chem.MolFromSmarts("[C]1[C]([O])[C][C][C]1")
        if not mol.HasSubstructMatch(hydroxyl_on_ring):
            # Check for a THF ring for thromboxanes
            thf_ring_pattern = Chem.MolFromSmarts("C1CCOC1")
            if not mol.HasSubstructMatch(thf_ring_pattern):
                return False, "Has cyclopentane ring but does not match prostaglandin/thromboxane pattern"

    # Check for epoxide pattern: C-O-C
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")
    if mol.HasSubstructMatch(epoxide_pattern):
         pass # We've already checked for double bonds
        
    #Check for hydroperoxide
    hydroperoxide_pattern = Chem.MolFromSmarts("[O]-[O]")
    if mol.HasSubstructMatch(hydroperoxide_pattern):
        pass
    
    # Check for Leukotriene structures.
    # Look for conjugated double bonds and an alcohol function or epoxy in a chain
    # leukotriene_pattern = Chem.MolFromSmarts("C=CC=CC[C](O)")
    # if mol.HasSubstructMatch(leukotriene_pattern):
    #      pass
    # A more specific pattern. Includes a thioether.
    #leukotriene_pattern = Chem.MolFromSmarts("[C][S][C][C](O)[C]=C[C]=C[C]=C[C]")
    #if mol.HasSubstructMatch(leukotriene_pattern):
    #     pass
    
    
    # Special Case handling: Some very specific structures like latanoprost or some of the diglycerides are hard to match generically.
    # if smiles == "O[C@H]1[C@@H]([C@H]([C@@H](O)C1)C/C=C\C(C(CC(O)=O)([2H])[2H])([2H])[2H])CC[C@@H](O)CCC2=CC=CC=C2": # latanoprost-d4
    #    return True, "Matches Latanoprost-d4 (special case)"
    # if smiles == "O[C@H](CC[C@H]1[C@H](O)C[C@H](O)[C@@H]1C\C=C/CCCC(O)=O)CCc1ccccc1": # latanoprost
    #    return True, "Matches Latanoprost (special case)"
    # if smiles == "S1C=C(OCC(O)\C=C\[C@@H]2[C@H]([C@@H](O)C[C@H]2O)C/C=C\CCCC(O)=O)C=C1": # Tiaprost
    #    return True, "Matches Tiaprost (special case)"

    return True, "Matches Icosanoid criteria (C20 backbone with double bonds and oxygenation)"