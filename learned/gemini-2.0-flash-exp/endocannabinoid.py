"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is likely an endocannabinoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely an endocannabinoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for long carbon chain (fatty acid or similar)
    fatty_acid_pattern = Chem.MolFromSmarts("C[C,C=C]~[C,C=C]~[C,C=C]~[C,C=C]~[C,C=C]") # allows for unsaturated and saturated chains, more generic than the previous one.
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No characteristic fatty acid chain found"

    fatty_chain_carbons = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C[C,C=C]~[C,C=C]~[C,C=C]~[C,C=C]~[C,C=C]"))[0]) # get number of carbons of longest chain
    if fatty_chain_carbons < 10:
      return False, "Chain too short to be a fatty acid"

    # 2. Check for polar head group (ethanolamide, glycerol or similar)
    ethanolamide_pattern = Chem.MolFromSmarts("NCCO")
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    ether_pattern = Chem.MolFromSmarts("COC")


    has_ethanolamide = mol.HasSubstructMatch(ethanolamide_pattern)
    has_glycerol = mol.HasSubstructMatch(glycerol_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_amide = mol.HasSubstructMatch(amide_pattern)
    has_ether = mol.HasSubstructMatch(ether_pattern)


    if not (has_ethanolamide or has_glycerol) and not (has_ester or has_amide or has_ether):
         return False, "No polar head group and linker found"

    #3. Check for connection between polar head and fatty acid chain
    head_group_found= False
    if has_amide:
      amide_matches = mol.GetSubstructMatches(amide_pattern)
      for match in amide_matches: #iterates over each tuple of indices for each match
         for atom_index in match: #iterates over each atom index within a match
           neighboring_carbons= [atom.GetAtomicNum() for atom in mol.GetAtomWithIdx(atom_index).GetNeighbors() if atom.GetAtomicNum() == 6 ]
           if len(neighboring_carbons) > 0: # if connected to at least a carbon
              head_group_found = True
              break
         if head_group_found:
            break
    if has_ester and not head_group_found:
      ester_matches = mol.GetSubstructMatches(ester_pattern)
      for match in ester_matches:
        for atom_index in match:
            neighboring_carbons = [atom.GetAtomicNum() for atom in mol.GetAtomWithIdx(atom_index).GetNeighbors() if atom.GetAtomicNum() == 6]
            if len(neighboring_carbons) > 0:
                head_group_found = True
                break
        if head_group_found:
            break

    if has_ether and not head_group_found:
        ether_matches = mol.GetSubstructMatches(ether_pattern)
        for match in ether_matches:
             for atom_index in match:
                 neighboring_carbons = [atom.GetAtomicNum() for atom in mol.GetAtomWithIdx(atom_index).GetNeighbors() if atom.GetAtomicNum() == 6]
                 if len(neighboring_carbons) > 0:
                     head_group_found = True
                     break
             if head_group_found:
                break

    if not head_group_found and (has_ester or has_amide or has_ether):
         return False, "No link to fatty chain found"

    #4. Check for presence of an epoxy group (common in some endocannabinoids)
    epoxy_pattern = Chem.MolFromSmarts("C1OC1")
    has_epoxy = mol.HasSubstructMatch(epoxy_pattern)
    if has_epoxy:
        epoxy_matches = mol.GetSubstructMatches(epoxy_pattern)
        epoxy_carbon_found = False
        for match in epoxy_matches:
            for atom_index in match:
                neighboring_carbons = [atom.GetAtomicNum() for atom in mol.GetAtomWithIdx(atom_index).GetNeighbors() if atom.GetAtomicNum() == 6 ]
                if len(neighboring_carbons) > 0 :
                    epoxy_carbon_found=True
                    break
            if epoxy_carbon_found:
                break
        if not epoxy_carbon_found:
            return False, "Epoxy ring is not part of the fatty chain"

    # If all conditions met, it is likely an endocannabinoid
    return True, "Likely an endocannabinoid based on long chain, polar head group, and linker, epoxy group detected"