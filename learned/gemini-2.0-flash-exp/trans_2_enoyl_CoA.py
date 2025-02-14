"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    A trans-2-enoyl-CoA is a coenzyme A molecule with a fatty acid chain
    attached via a thioester bond with a trans double bond at the 2-position
    relative to the carbonyl.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trans-2-enoyl-CoA, False otherwise.
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the CoA SMARTS pattern. This is the minimal core of the CoA moiety
    # we are checking to avoid issues with tautomers.
    coa_pattern = Chem.MolFromSmarts('SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12')
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Find all C=C double bonds
    double_bond_pattern = Chem.MolFromSmarts("[C]=[C]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)

    found_correct_bond = False
    for match in double_bond_matches:
        carbon1_idx = match[0]
        carbon2_idx = match[1]

        # Check if the double bond is trans
        
        chiral_centers = Chem.FindMolChiralCenters(mol,force=True, includeUnassigned=True)

        is_trans = False

        for center in chiral_centers:
            if center[0] == carbon1_idx:
                if center[1] == 'E':
                  is_trans = True
            if center[0] == carbon2_idx:
                if center[1] == 'E':
                    is_trans = True

        if not is_trans:
          continue
          
        
        # Check for correct position of the double bond, next to a carbonyl.
        carbon1 = mol.GetAtomWithIdx(carbon1_idx)
        carbon2 = mol.GetAtomWithIdx(carbon2_idx)
        
        for neighbor1 in carbon1.GetNeighbors():
          if neighbor1.GetHybridization() == Chem.HybridizationType.SP2:
            for neighbor2 in neighbor1.GetNeighbors():
                if neighbor2.GetAtomicNum() == 8 and neighbor2.GetBond().GetBondType() == Chem.BondType.DOUBLE:
                    for thio_neighbor in neighbor1.GetNeighbors():
                       if thio_neighbor.GetAtomicNum() == 16:
                          found_correct_bond = True
                          break
                    if found_correct_bond:
                       break
            if found_correct_bond:
                break
            
        for neighbor1 in carbon2.GetNeighbors():
          if neighbor1.GetHybridization() == Chem.HybridizationType.SP2:
            for neighbor2 in neighbor1.GetNeighbors():
                if neighbor2.GetAtomicNum() == 8 and neighbor2.GetBond().GetBondType() == Chem.BondType.DOUBLE:
                    for thio_neighbor in neighbor1.GetNeighbors():
                       if thio_neighbor.GetAtomicNum() == 16:
                          found_correct_bond = True
                          break
                    if found_correct_bond:
                       break
            if found_correct_bond:
                break


        if found_correct_bond:
           break


    if not found_correct_bond:
        return False, "Could not find correct trans double bond at position 2"
    
    # Check molecular weight - generally > 700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, "Molecular weight too low for trans-2-enoyl-CoA"

    # Check number of carbon atoms in fatty acid chain
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
       return False, "Too few carbons in the fatty acid chain"
    

    return True, "Molecule is a trans-2-enoyl-CoA"