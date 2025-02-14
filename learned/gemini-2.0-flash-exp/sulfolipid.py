"""
Classifies: CHEBI:61384 sulfolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid contains a sulfonic acid residue linked by a carbon-sulfur bond to a lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for the sulfonate group (-OS(=O)(=O)O-) or sulfonic acid (-S(=O)(=O)OH) with a C-S bond.
    # Note the use of [CX4] to specify a carbon, and ~ to specify a single bond.
    sulfonate_pattern = Chem.MolFromSmarts("[CX4]~[S](=[OX1])(=[OX1])[OX2]")
    sulfonic_acid_pattern = Chem.MolFromSmarts("[CX4]~[S](=[OX1])(=[OX1])[OX1H0]")

    sulfonate_matches = mol.GetSubstructMatches(sulfonate_pattern)
    sulfonic_acid_matches = mol.GetSubstructMatches(sulfonic_acid_pattern)
    if not sulfonate_matches and not sulfonic_acid_matches:
       return False, "No sulfonate or sulfonic acid group with a C-S bond found."


    # Get atoms of the carbons directly bonded to sulfur
    sulfur_atoms = []
    if sulfonate_matches:
        for match in sulfonate_matches:
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 16:  # Sulfur
                    sulfur_atoms.append(atom)
    elif sulfonic_acid_matches:
        for match in sulfonic_acid_matches:
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 16:  # Sulfur
                    sulfur_atoms.append(atom)


    # Verify C-S Bond and find the directly connected Carbon atoms
    connected_carbons = []
    for sulfur_atom in sulfur_atoms:
      for neighbor in sulfur_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:
            bond = mol.GetBondBetweenAtoms(sulfur_atom.GetIdx(), neighbor.GetIdx())
            if bond and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
              connected_carbons.append(neighbor)
    if not connected_carbons:
      return False, "No carbon directly bonded to sulfur via a single bond."



    # 3. Check if the carbon is part of a lipid-like moiety
    is_lipid = False
    for carbon_atom in connected_carbons:
        # Check for long fatty acid chains with at least 8 carbons
        fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
        match_list = mol.GetSubstructMatches(fatty_acid_pattern)
        for match in match_list:
          if carbon_atom.GetIdx() in match:
            is_lipid = True
            break
        if is_lipid:
            break

        # Check for a fatty acid connected through an ester or amide
        acyl_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3](=[OX1])[OX2,NX2]")
        match_list = mol.GetSubstructMatches(acyl_pattern)
        for match in match_list:
          if carbon_atom.GetIdx() in match:
            is_lipid = True
            break
        if is_lipid:
            break
        # Check for common sugar residues
        sugar_pattern = Chem.MolFromSmarts("[OX2][CX4][CX4][CX4][CX4][CX4][OX2]") #simplified sugar ring
        match_list = mol.GetSubstructMatches(sugar_pattern)
        for match in match_list:
            if carbon_atom.GetIdx() in match:
                is_lipid = True
                break
        if is_lipid:
            break

        # Check for glycosidic linkages
        glycosidic_pattern = Chem.MolFromSmarts("[OX2][CX4][CX4][CX4][CX4][CX4][OX2]~[CX4,CX3]")
        match_list = mol.GetSubstructMatches(glycosidic_pattern)
        for match in match_list:
            if carbon_atom.GetIdx() in match:
                is_lipid = True
                break
        if is_lipid:
          break

    if not is_lipid:
      return False, "Carbon connected to sulfur is not part of a recognizable lipid-like moiety"


    # 4. Check molecular weight to discard obvious non-lipids
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 450:
      return False, "Molecular weight too low for a sulfolipid"



    return True, "Contains a sulfonic acid residue linked to a lipid via a carbon-sulfur bond."