"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    Waxes are generally long-chain esters or alcohols with flexibility
    and moderate to high molecular weight.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 0. Check for heteroatoms
    heteroatoms = [atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6, 8]]
    if len(heteroatoms) > 0:
        return False, "Contains heteroatoms other than C, H and O; not typical for wax"

    # 1. Check for long continuous carbon chains (at least 16) - more specific check
    
    # Get all carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]

    longest_chain = 0
    for start_atom in carbon_atoms:
        for end_atom in carbon_atoms:
            if start_atom == end_atom:
              continue
            
            path = Chem.GetShortestPath(mol, start_atom.GetIdx(), end_atom.GetIdx())
            
            if path is not None and len(path)>longest_chain: # path exists, and longer than longest
                is_all_carbon = all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in path)
                if is_all_carbon:
                    longest_chain = len(path)
    
    if longest_chain < 15:
        return False, f"Too few carbons in the longest chain ({longest_chain}); typical waxes should have at least 15 contiguous C atoms"
    
    # Count total number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

     # 2. Check for ester linkages (or long-chain alcohols)
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4]") # Ester next to a carbon
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H1]") # long-chain alcohols check
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)

    if len(ester_matches) == 0 and len(alcohol_matches) == 0:
        return False, f"No long-chain ester or alcohol groups detected."

    # Check that the ester/alcohol group is part of the long chain. If it's an ester it must have another alkyl
    if len(ester_matches) > 0 :
        ester_match_ok = False
        for match in ester_matches:
           
           ester_carbon1 = mol.GetAtomWithIdx(match[0])
           ester_carbon2 = mol.GetAtomWithIdx(match[3]) # carbon to which the ester oxygen is linked
           if (ester_carbon1 in carbon_atoms and ester_carbon2 in carbon_atoms):
              ester_match_ok = True
              break
        if not ester_match_ok:
             return False, "Ester group is not part of long carbon chain"
            
    if len(alcohol_matches) > 0:
        alcohol_match_ok = False
        for match in alcohol_matches:
            alcohol_carbon = mol.GetAtomWithIdx(match[0])
            if (alcohol_carbon in carbon_atoms):
                alcohol_match_ok = True
                break
        if not alcohol_match_ok:
             return False, "Alcohol group is not part of long carbon chain"

     #3 Correct check number of oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    required_o_for_esters = 0
    if len(ester_matches) > 0:
       required_o_for_esters = len(ester_matches)*2
    required_o_for_alcohols = 0
    if len(alcohol_matches) > 0:
      required_o_for_alcohols = len(alcohol_matches)

    if o_count != required_o_for_esters + required_o_for_alcohols:
         return False, f"Inconsistent number of oxygens {o_count}, expected {required_o_for_esters+required_o_for_alcohols} based on ester/alcohol count."
   

    # 4. Check for rotatable bonds (at least 8)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, f"Too few rotatable bonds ({n_rotatable}). Should be at least 8"

    # 5. Check molecular weight (approx range 300-1000)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1000:
        return False, f"Molecular weight out of range ({mol_wt}). Should be between 300 and 1000"

    # 6. Check the carbon ratio
    chain_ratio = longest_chain / c_count
    if chain_ratio < 0.7:
        return False, "Ratio of carbons in longest chain to total carbons is too low"
    
    return True, "Meets the criteria for a wax (long, largely saturated chains, ester or alcohol, flexible, good MW)."