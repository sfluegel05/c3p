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

    # Check for carboxylic acid
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group"

    # Check for long chain (16-20 carbons backbone)
    long_chain_pattern = Chem.MolFromSmarts("C[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Not a long-chain fatty acid"
    
    # Check for prostanoid ring system
    prostanoid_pattern = Chem.MolFromSmarts("C1CC[C]2[C]1[C](C[C]([C]2)O)[C]=O") #simplified for a quick check
    if mol.HasSubstructMatch(prostanoid_pattern):
        return False, "Likely a Prostanoid (contains characteristic ring)"


    # Check for conjugated triene pattern (leukotriene-like) with more flexibility for branching
    # this pattern now captures C=C-C=C-C=C and C=C-C=C-C=C-C with some flexibility
    leukotriene_pattern = Chem.MolFromSmarts("[C]=[C]-[C]=[C]-[C]=[C]") # basic triene

    if mol.HasSubstructMatch(leukotriene_pattern):
           # further check to discard molecules where the triene part is not part of the main chain
           if mol.GetSubstructMatch(Chem.MolFromSmarts("[C]=[C]-[C]=[C]-[C]=[C]-[C]")):
              return False, "Likely a Leukotriene (conjugated triene pattern)"

    # Check for lipoxin/resolvin features - multiple double bonds and OH/epoxide/carbonyl
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    hydroxyl_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("[OH1]"))
    epoxide_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("C1OC1"))
    carbonyl_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("C=O"))

    # Oxygenation pattern check
    if (double_bond_count >= 2 and (len(hydroxyl_groups) >= 2 or len(epoxide_groups) >= 1)):
        pass #likely a non-classic icosanoid
    elif (double_bond_count >= 3 and len(hydroxyl_groups) >= 1):
        pass
    else:
      return False, "Not enough double bonds and hydroxyl/epoxide groups"

    # Molecular Weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low"
    
    return True, "Likely a non-classic icosanoid"