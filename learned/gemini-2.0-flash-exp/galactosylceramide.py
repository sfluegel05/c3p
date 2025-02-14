"""
Classifies: CHEBI:36498 galactosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide consists of a galactose head group linked to a ceramide.
    Ceramides have a sphingosine or sphinganine backbone with an amide-linked fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a galactosylceramide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for alpha and beta galactopyranose (C1 position specified)
    alpha_galactose_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)O)")
    beta_galactose_pattern = Chem.MolFromSmarts("[C@@H]1([C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)O)")
    
    #Check for a galactopyranose moiety.
    if not (mol.HasSubstructMatch(alpha_galactose_pattern) or mol.HasSubstructMatch(beta_galactose_pattern)):
            return False, "No alpha or beta galactose moiety found"
    
    # SMARTS for sphingosine/sphinganine backbone including double bond variation and hydroxyl groups
    # Includes the characteristic 2-amino-1,3-diol
    sphingosine_core = Chem.MolFromSmarts("[C@H]([C@H](O)[C@H]([C])O)NC(=O)")  
    phytosphingosine_core = Chem.MolFromSmarts("[C@H]([C@H](O)[C@H]([C])(O)O)NC(=O)")

    if not (mol.HasSubstructMatch(sphingosine_core) or mol.HasSubstructMatch(phytosphingosine_core)):
        return False, "No sphingosine/phytosphingosine backbone found"

    # Targeted glycosidic bond check: Oxygen between galactose C1 and ceramide backbone carbon
    # This enforces connectivity between the two, ensuring a galactosyl-ceramide bond is detected.
    
    galactose_c1_pattern = Chem.MolFromSmarts("[C:1][OH1][C@H]([C@H]([C@H]([C@H](O[*:2])CO)O)O)O")
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C:2][OX2][C:3]") # C:2 must be ceramide carbon
    
    galactose_matches = mol.GetSubstructMatches(galactose_c1_pattern)
    if not galactose_matches:
        return False, "No glycosidic bond present between galactose and ceramide"

    glycosidic_bond_present = False
    for gal_match in galactose_matches:
        gal_c1 = gal_match[0]
        for atom in mol.GetAtomWithIdx(gal_c1).GetNeighbors():
            if atom.GetAtomicNum() == 8: # Oxygen
                for ceramide_neighbor in atom.GetNeighbors():
                        if ceramide_neighbor.HasSubstructMatch(sphingosine_core) or ceramide_neighbor.HasSubstructMatch(phytosphingosine_core):
                            glycosidic_bond_present = True
                            break
        if glycosidic_bond_present:
            break
            
    if not glycosidic_bond_present:
      return False, "No glycosidic bond present between galactose and ceramide"
    
    # Fatty acid chain check
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "Missing fatty acid chain"
    

    return True, "Contains a galactose moiety linked to a ceramide backbone."