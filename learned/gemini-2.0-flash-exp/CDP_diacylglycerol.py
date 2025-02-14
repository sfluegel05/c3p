"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    A CDP-diacylglycerol has a cytidine diphosphate group linked to a diacylglycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a CDP-diacylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CDP core.
    # Core pattern: Cytosine-Ribose-Diphosphate
    cdp_core_pattern = Chem.MolFromSmarts('C1[C@H]([C@H]([C@@H]([C@H](O1)O)N2C=CC(=NC2=O)N)O)COP(=O)(O)OP(=O)(O)O')
    if not mol.HasSubstructMatch(cdp_core_pattern):
       return False, "CDP core not found"
    
    # Look for the glycerol phosphate pattern connected to the CDP core:
    # This pattern identifies glycerol attached to diphosphate.
    glycerol_phosphate_pattern = Chem.MolFromSmarts("OC[C@H](CO[P](=O)(O)O[P](=O)(O)O)[CH2]O")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "Glycerol-phosphate part of CDP-DG not found"

    # Look for 2 ester groups connected to the glycerol
    ester_pattern = Chem.MolFromSmarts("[CX4]([OX2])[CX3](=[OX1])~[CX4]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    #  Count ester groups in the right position - i.e. adjacent to the glycerol.
    
    count = 0
    for match in ester_matches:
            atom_idx = match[0] # The oxygen atom
            
            #Check if the oxygen is attached to the glycerol part
            
            for atom in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                if atom.GetSymbol() == "C" and atom.GetDegree() == 4:
                    glycerol_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("OC[C@H](CO[P](=O)(O)O[P](=O)(O)O)[CH2]O"))
                    for gly_match in glycerol_matches:
                        for gly_atom_idx in gly_match:
                           if mol.GetAtomWithIdx(gly_atom_idx).GetSymbol()=="C" and mol.GetAtomWithIdx(gly_atom_idx).GetDegree() == 4:
                                   if atom.GetIdx() == mol.GetAtomWithIdx(gly_atom_idx).GetIdx():
                                     count += 1
                                     break
                        if count > 2:
                            break    
            if count > 2:
              break
    if count < 2:
        return False, f"Missing at least two acyl groups, got {count}"
        

    # Count phosphorus and nitrogen to ensure correct stoichiometry.
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if p_count != 2:
        return False, "Must have exactly 2 phosphorus atoms (diphosphate group)"
    if n_count != 3:
         return False, "Must have 3 nitrogen atoms in CDP base"

    return True, "Contains CDP core with diacylglycerol attached"