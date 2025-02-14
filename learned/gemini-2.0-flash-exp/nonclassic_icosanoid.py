"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a non-classic icosanoid based on its SMILES string.
    A non-classic icosanoid is an oxygenated C20 fatty acid that is NOT a leukotriene or prostanoid.
    Lipoxins and resolvins are included.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple(bool, str): (True, reason) if molecule is a non-classic icosanoid, (False, reason) otherwise
                         or (None, None) if the SMILES string is invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None

    # 1. Check for C20 fatty acid backbone with carboxyl group
    # This SMARTS pattern checks for 20 carbons in a chain and a terminal carboxylic acid
    c20_acid_pattern = Chem.MolFromSmarts("[CX4]([CX4])([CX4])([CX4])~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(c20_acid_pattern):
      return False, "Not a C20 fatty acid"
    
    # 2. Check number of carbons to confirm C20 origin
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
       return False, "Does not contain 20 carbons"
   
    # 3. Check for prostanoid ring system
    prostanoid_pattern = Chem.MolFromSmarts("C1CC[C]2[C]1[C](C[C]([C]2)O)[C]=O") #simplified for a quick check
    if mol.HasSubstructMatch(prostanoid_pattern):
        return False, "Likely a Prostanoid (contains characteristic ring)"

    # 4. Check for leukotriene conjugated triene pattern in main chain
    #    Leukotrienes are usually linear and have three conjugated double bonds within the main chain.
    #    This check is now more permissive of branching, but needs to be in the main chain.
    main_chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX3]")
    main_chain_match = mol.GetSubstructMatch(main_chain_pattern)
    if main_chain_match:
      main_chain = Chem.PathToSubmol(mol,main_chain_match)
      leukotriene_pattern = Chem.MolFromSmarts("[C]=[C]-[C]=[C]-[C]=[C]") #triene pattern
      if main_chain.HasSubstructMatch(leukotriene_pattern):
        return False, "Likely a Leukotriene (conjugated triene pattern in main chain)"

      conjugated_tetraene = Chem.MolFromSmarts("[C]=[C]-[C]=[C]-[C]=[C]-[C]=[C]") # tetraene pattern
      if main_chain.HasSubstructMatch(conjugated_tetraene):
        return False, "Likely a Leukotriene (conjugated tetraene pattern in main chain)"

    # 5. Check for oxygenation features (double bonds, hydroxyls, epoxides, carbonyls)
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    hydroxyl_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("[OH1]"))
    epoxide_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("C1OC1"))
    carbonyl_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("C=O"))

    # Oxygenation pattern check (simplified from the previous version as C20 check is now more strict)
    if (double_bond_count >= 2 and (len(hydroxyl_groups) >= 2 or len(epoxide_groups) >= 1)):
        pass #likely a non-classic icosanoid
    elif (double_bond_count >= 3 and len(hydroxyl_groups) >= 1):
        pass
    else:
      return False, "Not enough double bonds and hydroxyl/epoxide groups"

    # 6. Molecular Weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low"

    # 7. Count oxygens to verify oxygenation pattern
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Too few oxygens, should be at least 2"
    
    return True, "Likely a non-classic icosanoid"